// 模块功能：实现二维砌体宏单元的非核心框架，包括命令解析、节点连接、状态管理、响应和序列化。
// 使用流程：先定义三类单轴材料和二维 geomTransf，再调用 element masonryMacro 创建本单元。
// 输入输出：输入节点、材料和坐标变换；输出经过内部剪切变形凝聚后的刚度与恢复力。

#include "MasonryMacro2D.h"
#include "MasonryBendingMat.h"
#include "MasonryShearMat.h"

#include <Channel.h>
#include <CrdTransf.h>
#include <Domain.h>
#include <ElementResponse.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <Renderer.h>
#include <UniaxialMaterial.h>
#include <elementAPI.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace {

// material: 需要分配数据库标签的材料。
// theChannel: 提供数据库标签的通讯通道。
int ensureMaterialDbTag(UniaxialMaterial *material, Channel &theChannel)
{
    int dbTag = material->getDbTag();
    if (dbTag == 0) {
        dbTag = theChannel.getDbTag();
        material->setDbTag(dbTag);
    }
    return dbTag;
}

// material: 需要接收或重建的材料指针。
// classTag: 目标材料类型标签。
// dbTag: 目标材料数据库标签。
// commitTag: 当前提交标签。
// theChannel: 提供材料数据的通讯通道。
// theBroker: 创建材料对象的工厂。
int receiveMaterial(UniaxialMaterial *&material, int classTag, int dbTag, int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (material == 0 || material->getClassTag() != classTag) {
        delete material;
        material = theBroker.getNewUniaxialMaterial(classTag);
    }
    if (material == 0) {
        return -1;
    }
    material->setDbTag(dbTag);
    return material->recvSelf(commitTag, theChannel, theBroker);
}

// material: 采用弹性路径试算、且不推进非线性历史的砌体材料。
// strain: 当前试算变形。
int setTrialLinearStrain(UniaxialMaterial *material, double strain)
{
    if (MasonryBendingMat *bendingMaterial = dynamic_cast<MasonryBendingMat *>(material)) {
        return bendingMaterial->setTrialLinearStrain(strain);
    }
    if (MasonryShearMat *shearMaterial = dynamic_cast<MasonryShearMat *>(material)) {
        return shearMaterial->setTrialLinearStrain(strain);
    }
    return -1;
}

// material: 需要读取当前骨架屈服变形的砌体材料。
double getYieldStrain(UniaxialMaterial *material)
{
    if (MasonryBendingMat *bendingMaterial = dynamic_cast<MasonryBendingMat *>(material)) {
        return bendingMaterial->getYieldStrain();
    }
    if (MasonryShearMat *shearMaterial = dynamic_cast<MasonryShearMat *>(material)) {
        return shearMaterial->getYieldStrain();
    }
    return 0.0;
}

} // namespace

// 解析二维 masonryMacro 命令并创建单元。
// 命令格式：element masonryMacro tag iNode jNode bendingMat shearMat axialMat transfTag
//           <-AxialForceInteraction width thickness fm cohesion ft <kd> <mu>>
//           <-ExclusiveFailureMode>
//           <-LocalIteration maximumIterations relativeTolerance>
//           <-LocalIterationDisplacement maximumIterations relativeTolerance>
void *OPS_MasonryMacro2D(void)
{
    const int numberOfArguments = OPS_GetNumRemainingInputArgs();
    if (numberOfArguments < 7) {
        opserr << "WARNING incorrect arguments: element masonryMacro tag iNode jNode bendingMat shearMat axialMat transfTag <-ExclusiveFailureMode> <-AxialForceInteraction width thickness fm cohesion ft <kd> <mu>> <-LocalIteration maximumIterations relativeTolerance> <-LocalIterationDisplacement maximumIterations relativeTolerance>" << endln;
        return 0;
    }

    int integerData[7];
    int numberOfData = 7;
    if (OPS_GetIntInput(&numberOfData, integerData) < 0) {
        opserr << "WARNING masonryMacro2D failed to read integer arguments" << endln;
        return 0;
    }

    UniaxialMaterial *bendingMaterial = OPS_getUniaxialMaterial(integerData[3]);
    UniaxialMaterial *shearMaterial = OPS_getUniaxialMaterial(integerData[4]);
    UniaxialMaterial *axialMaterial = OPS_getUniaxialMaterial(integerData[5]);
    CrdTransf *coordinateTransformation = OPS_getCrdTransf(integerData[6]);

    if (bendingMaterial == 0 || shearMaterial == 0 || axialMaterial == 0 || coordinateTransformation == 0) {
        opserr << "WARNING masonryMacro2D could not find a material or coordinate transformation for element " << integerData[0] << endln;
        return 0;
    }

    // 未输入截面和材料参数时保持原行为，不启用轴力相关骨架更新。
    double interactionData[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.85, 0.4};
    int maximumIterations = 50;
    double relativeTolerance = 1.0e-8;
    bool useDisplacementConvergence = true;
    bool useExclusiveFailureMode = false;
    bool localIterationOptionSpecified = false;
    int remainingArguments = numberOfArguments - 7;
    while (remainingArguments > 0) {
        const char *option = OPS_GetString();
        --remainingArguments;

        if (std::strcmp(option, "-ExclusiveFailureMode") == 0) {
            useExclusiveFailureMode = true;
            continue;
        }

        if (std::strcmp(option, "-LocalIteration") == 0 || std::strcmp(option, "-LocalIterationDisplacement") == 0) {
            if (localIterationOptionSpecified) {
                opserr << "WARNING masonryMacro2D -LocalIteration and -LocalIterationDisplacement are mutually exclusive" << endln;
                return 0;
            }
            const bool displacementOption = std::strcmp(option, "-LocalIterationDisplacement") == 0;
            if (remainingArguments < 2) {
                opserr << "WARNING masonryMacro2D " << option << " requires maximumIterations and relativeTolerance" << endln;
                return 0;
            }
            int numberOfData = 1;
            if (OPS_GetIntInput(&numberOfData, &maximumIterations) < 0 || OPS_GetDoubleInput(&numberOfData, &relativeTolerance) < 0) {
                opserr << "WARNING masonryMacro2D failed to read " << option << " parameters" << endln;
                return 0;
            }
            remainingArguments -= 2;
            useDisplacementConvergence = displacementOption;
            localIterationOptionSpecified = true;
            continue;
        }

        if (std::strcmp(option, "-AxialForceInteraction") != 0) {
            opserr << "WARNING masonryMacro2D unknown option " << option << endln;
            return 0;
        }

        // 若后面还有局部迭代选项，则为它保留三个输入；否则剩余输入全部属于轴力相关参数。
        int numberOfInteractionData = remainingArguments > 7 ? remainingArguments - 3 : remainingArguments;
        if (numberOfInteractionData < 5 || numberOfInteractionData > 7) {
            opserr << "WARNING masonryMacro2D -AxialForceInteraction requires width thickness fm cohesion ft <kd> <mu>" << endln;
            return 0;
        }
        if (OPS_GetDoubleInput(&numberOfInteractionData, interactionData) < 0) {
            opserr << "WARNING masonryMacro2D failed to read axial-interaction parameters" << endln;
            return 0;
        }
        remainingArguments -= numberOfInteractionData;
    }

    return new MasonryMacro2D(integerData[0], integerData[1], integerData[2], *bendingMaterial, *shearMaterial, *axialMaterial, *coordinateTransformation, interactionData[0], interactionData[1], interactionData[2], interactionData[3], interactionData[4], interactionData[5], interactionData[6], maximumIterations, relativeTolerance, useDisplacementConvergence, useExclusiveFailureMode);
}

MasonryMacro2D::MasonryMacro2D(int tag, int nodeI, int nodeJ, UniaxialMaterial &bendingMaterial, UniaxialMaterial &shearMaterialInput, UniaxialMaterial &axialMaterialInput, CrdTransf &coordinateTransformation, double widthInput, double thicknessInput, double compressiveStrengthInput, double cohesionInput, double diagonalTensileStrengthInput, double stressBlockCoefficientInput, double frictionCoefficientInput, int maximumIterationsInput, double relativeToleranceInput, bool useDisplacementConvergenceInput, bool useExclusiveFailureModeInput)
    : Element(tag, ELE_TAG_MasonryMacro2D), connectedExternalNodes(2), theCoordTransf(coordinateTransformation.getCopy2d()), shearMaterial(shearMaterialInput.getCopy()), axialMaterial(axialMaterialInput.getCopy()), width(widthInput), thickness(thicknessInput), compressiveStrength(compressiveStrengthInput), cohesion(cohesionInput), diagonalTensileStrength(diagonalTensileStrengthInput), stressBlockCoefficient(stressBlockCoefficientInput), frictionCoefficient(frictionCoefficientInput), committedContraflexureDistance(0.0), maximumIterations(maximumIterationsInput), relativeTolerance(relativeToleranceInput), useDisplacementConvergence(useDisplacementConvergenceInput), useExclusiveFailureMode(useExclusiveFailureModeInput), committedFailureMode(0), trialFailureMode(0), tangentStiffness(6, 6), initialStiffness(6, 6), resistingForce(6)
{
    connectedExternalNodes(0) = nodeI;
    connectedExternalNodes(1) = nodeJ;
    theNodes[0] = 0;
    theNodes[1] = 0;
    bendingMaterials[0] = bendingMaterial.getCopy();
    bendingMaterials[1] = bendingMaterial.getCopy();
}

MasonryMacro2D::MasonryMacro2D()
    : Element(0, ELE_TAG_MasonryMacro2D), connectedExternalNodes(2), theCoordTransf(0), shearMaterial(0), axialMaterial(0), width(0.0), thickness(0.0), compressiveStrength(0.0), cohesion(0.0), diagonalTensileStrength(0.0), stressBlockCoefficient(0.85), frictionCoefficient(0.4), committedContraflexureDistance(0.0), maximumIterations(50), relativeTolerance(1.0e-8), useDisplacementConvergence(true), useExclusiveFailureMode(false), committedFailureMode(0), trialFailureMode(0), tangentStiffness(6, 6), initialStiffness(6, 6), resistingForce(6)
{
    theNodes[0] = 0;
    theNodes[1] = 0;
    bendingMaterials[0] = 0;
    bendingMaterials[1] = 0;
}

MasonryMacro2D::~MasonryMacro2D()
{
    delete theCoordTransf;
    delete bendingMaterials[0];
    delete bendingMaterials[1];
    delete shearMaterial;
    delete axialMaterial;
}

int MasonryMacro2D::getNumExternalNodes(void) const
{
    return 2;
}

const ID &MasonryMacro2D::getExternalNodes(void)
{
    return connectedExternalNodes;
}

Node **MasonryMacro2D::getNodePtrs(void)
{
    return theNodes;
}

int MasonryMacro2D::getNumDOF(void)
{
    return 6;
}

void MasonryMacro2D::setDomain(Domain *theDomain)
{
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        this->DomainComponent::setDomain(0);
        return;
    }

    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));
    if (theNodes[0] == 0 || theNodes[1] == 0) {
        opserr << "MasonryMacro2D::setDomain - element " << this->getTag() << " could not find both nodes" << endln;
        return;
    }
    if (theNodes[0]->getNumberDOF() != 3 || theNodes[1]->getNumberDOF() != 3) {
        opserr << "MasonryMacro2D::setDomain - element " << this->getTag() << " requires three DOFs at each node" << endln;
        return;
    }
    if (theCoordTransf == 0 || theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
        opserr << "MasonryMacro2D::setDomain - element " << this->getTag() << " failed to initialize coordinate transformation" << endln;
        return;
    }

    if (theCoordTransf->getInitialLength() <= 0.0) {
        opserr << "MasonryMacro2D::setDomain - element " << this->getTag() << " has zero length" << endln;
        return;
    }
    if (committedContraflexureDistance == 0.0) {
        committedContraflexureDistance = theCoordTransf->getInitialLength() / 2.0;
    }
    this->DomainComponent::setDomain(theDomain);
}

int MasonryMacro2D::commitState(void)
{
    const double momentI = bendingMaterials[0]->getStress();
    const double momentJ = bendingMaterials[1]->getStress();
    const double momentSum = momentI + momentJ;
    const double V = std::abs(momentSum / theCoordTransf->getInitialLength());
    double trialContraflexureDistance = committedContraflexureDistance;
    if (V > 0.0) {
        trialContraflexureDistance = std::min(std::max(std::abs(momentI), std::abs(momentJ)) / V, 1.0e8);
    }
    int result = this->Element::commitState();
    result += theCoordTransf->commitState();
    result += bendingMaterials[0]->commitState();
    result += bendingMaterials[1]->commitState();
    result += shearMaterial->commitState();
    result += axialMaterial->commitState();
    if (result == 0) {
        committedContraflexureDistance = trialContraflexureDistance;
        committedFailureMode = trialFailureMode;
    }
    return result;
}

int MasonryMacro2D::revertToLastCommit(void)
{
    int result = theCoordTransf->revertToLastCommit();
    result += bendingMaterials[0]->revertToLastCommit();
    result += bendingMaterials[1]->revertToLastCommit();
    result += shearMaterial->revertToLastCommit();
    result += axialMaterial->revertToLastCommit();
    trialFailureMode = committedFailureMode;
    return result;
}

int MasonryMacro2D::revertToStart(void)
{
    int result = theCoordTransf->revertToStart();
    result += bendingMaterials[0]->revertToStart();
    result += bendingMaterials[1]->revertToStart();
    result += shearMaterial->revertToStart();
    result += axialMaterial->revertToStart();
    committedContraflexureDistance = theCoordTransf->getInitialLength() / 2.0;
    committedFailureMode = 0;
    trialFailureMode = 0;
    resistingForce.Zero();
    return result;
}

int MasonryMacro2D::update(void)
{
    int result = theCoordTransf->update();
    if (result != 0) {
        return result;
    }

    // 从坐标变换读取节点基本位移，计算轴向、剪切和两端弯曲变形；然后将这些变形传递给对应的材料。
    //  I -------- J
    // 
    // 假设两端相对于轴线法线的转角为 thetaI 和 thetaJ （逆时针为正）
    // 两端节点弹簧的转角为 phiI 和 phiJ （逆时针为正），剪切弹簧的变形为 Delta （J节点在上为正）
    // 两端节点弯矩为 Mi 和 Mj （逆时针为正），两端节点剪力为 V （向上为正，垂直于刚性杆而不是轴线）
    // 轴力为 N （拉伸为正，注意不是轴向力，而是平行于刚性杆）；
    // 严格的轴向力为 N_ax = N * cos(alpha) + V * sin(alpha)，其中 alpha 为刚性杆与轴线的夹角。
    // sin(alpha) = Delta / L，L 为长度。但是这里近似认为 N_ax = N，忽略了剪切变形对轴向力的影响。
    // 
    // 那么可以列方程组（其中 phiI, phiJ, Delta 未知）：
    // 
    //  thetaI = phiI - Delta / L
    //  thetaJ = phiJ - Delta / L
    //  Mi(phiI) + Mj(phiJ) + V(Delta) * L =  N * Delta    (平衡方程)
    // 
    // 将上两个方程代入第三个方程，得到 (仅关于 Delta)：
    // Mi(thetaI + Delta / L) + Mj(thetaJ + Delta / L) + V(Delta) * L =  N * Delta

    // 获取当前单元的基本位移向量 [轴向，转角I，转角J]
    const Vector &eleTrialDeformation = theCoordTransf->getBasicTrialDisp();

    // 先更新轴向材料，获得当前平衡方程所需的轴力。
    const double axialDeformation = eleTrialDeformation(0);
    result = axialMaterial->setTrialStrain(axialDeformation);
    if (result != 0) {
        return result;
    }

    // 反弯点在本步固定，轴力仍采用当前试算值更新骨架。
    result = this->updateMaterialBackbones(axialMaterial->getStress());
    if (result != 0) {
        return result;
    }

    trialFailureMode = committedFailureMode;
    result = this->solveInternalShearDeformation(eleTrialDeformation(1), eleTrialDeformation(2), axialMaterial->getStress());
    if (result != 0 || !useExclusiveFailureMode || committedFailureMode != 0) {
        return result;
    }

    // 首次越过屈服点后锁定控制模式，并按锁定后的弹簧组合重新满足内部平衡。
    result = this->selectTrialFailureMode();
    if (result != 0 || trialFailureMode == 0) {
        return result;
    }
    return this->solveInternalShearDeformation(eleTrialDeformation(1), eleTrialDeformation(2), axialMaterial->getStress());
}

int MasonryMacro2D::updateMaterialBackbones(double axialForce)
{
    // 旧命令未提供相互作用参数时，不改变材料创建时输入的骨架。
    if (width == 0.0) {
        return 0;
    }

    MasonryBendingMat *bendingMaterialI = dynamic_cast<MasonryBendingMat *>(bendingMaterials[0]);
    MasonryBendingMat *bendingMaterialJ = dynamic_cast<MasonryBendingMat *>(bendingMaterials[1]);
    MasonryShearMat *masonryShearMaterial = dynamic_cast<MasonryShearMat *>(shearMaterial);
    if (bendingMaterialI == 0 || bendingMaterialJ == 0 || masonryShearMaterial == 0) {
        opserr << "MasonryMacro2D::updateMaterialBackbones - element " << this->getTag() << " requires MasonryBendingMat and MasonryShearMat when axial interaction is enabled" << endln;
        return -1;
    }

    // 单元内部轴力以拉伸为正；论文承载力公式中的 P 为正的轴向压力。
    const double compressionForce = std::max(0.0, -axialForce);

    // 论文式（7）：rocking/crushing 控制的弯矩承载力。
    const double compressionCapacity = stressBlockCoefficient * width * thickness * compressiveStrength;
    const double momentCapacity = compressionForce * width / 2.0 * (1.0 - compressionForce / compressionCapacity);

    // 论文式（13）：沿水平灰缝滑移的剪切承载力。
    double slidingCapacity = 0.0;
    if (compressionForce > 0.0) {
        slidingCapacity = (1.5 * width * thickness * cohesion + frictionCoefficient * compressionForce) / (1.0 + 3.0 * committedContraflexureDistance * thickness * cohesion / compressionForce);
    }

    // 论文式（14）：对角开裂剪切承载力，形状系数 xi 限制在 1.0～1.5。
    const double shapeFactor = std::min(std::max(theCoordTransf->getInitialLength() / width, 1.0), 1.5);
    const double tensileAreaForce = diagonalTensileStrength * width * thickness;
    const double diagonalCapacity = tensileAreaForce / shapeFactor * std::sqrt(1.0 + compressionForce / tensileAreaForce);
    const double shearCapacity = std::min(slidingCapacity, diagonalCapacity);

    int result = bendingMaterialI->setTrialBackbone(momentCapacity);
    result += bendingMaterialJ->setTrialBackbone(momentCapacity);
    result += masonryShearMaterial->setTrialBackbone(shearCapacity);
    return result;
}

int MasonryMacro2D::solveInternalShearDeformation(double thetaI, double thetaJ, double axialForce)
{
    const double initialShearDeformation = shearMaterial->getStrain();
    if (this->solveInternalShearDeformationByNewton(thetaI, thetaJ, axialForce, initialShearDeformation) == 0) {
        return 0;
    }
    if (this->solveInternalShearDeformationByBisection(thetaI, thetaJ, axialForce, initialShearDeformation) == 0) {
        return 0;
    }

    opserr << "MasonryMacro2D::solveInternalShearDeformation - element " << this->getTag() << " did not converge" << endln;
    return -1;
}

int MasonryMacro2D::evaluateInternalShearResidual(double thetaI, double thetaJ, double axialForce, double shearDeformation, double &residual, double &residualScale)
{
    const double length = theCoordTransf->getInitialLength();
    const double rotationI = thetaI + shearDeformation / length;
    const double rotationJ = thetaJ + shearDeformation / length;
    int result = 0;
    if (useExclusiveFailureMode && trialFailureMode == 0) {
        // 锁定前两类弹簧均按弹性直线试算，以弹性需求比判断真正首先屈服的模式。
        result = setTrialLinearStrain(shearMaterial, shearDeformation);
        result += setTrialLinearStrain(bendingMaterials[0], rotationI);
        result += setTrialLinearStrain(bendingMaterials[1], rotationJ);
    } else if (useExclusiveFailureMode && trialFailureMode == 1) {
        result = shearMaterial->setTrialStrain(shearDeformation);
        result += setTrialLinearStrain(bendingMaterials[0], rotationI);
        result += setTrialLinearStrain(bendingMaterials[1], rotationJ);
    } else if (useExclusiveFailureMode && trialFailureMode == 2) {
        result = setTrialLinearStrain(shearMaterial, shearDeformation);
        result += bendingMaterials[0]->setTrialStrain(rotationI);
        result += bendingMaterials[1]->setTrialStrain(rotationJ);
    } else {
        result = shearMaterial->setTrialStrain(shearDeformation);
        result += bendingMaterials[0]->setTrialStrain(rotationI);
        result += bendingMaterials[1]->setTrialStrain(rotationJ);
    }
    if (result != 0) {
        return result;
    }

    const double momentI = bendingMaterials[0]->getStress();
    const double momentJ = bendingMaterials[1]->getStress();
    const double shearForce = shearMaterial->getStress();
    residual = momentI + momentJ + shearForce * length - axialForce * shearDeformation;
    residualScale = std::max(1.0, std::abs(momentI) + std::abs(momentJ) + std::abs(shearForce * length) + std::abs(axialForce * shearDeformation));
    return 0;
}

int MasonryMacro2D::selectTrialFailureMode(void)
{
    const double shearYieldStrain = getYieldStrain(shearMaterial);
    const double bendingYieldStrainI = getYieldStrain(bendingMaterials[0]);
    const double bendingYieldStrainJ = getYieldStrain(bendingMaterials[1]);
    if (shearYieldStrain <= 0.0 || bendingYieldStrainI <= 0.0 || bendingYieldStrainJ <= 0.0) {
        opserr << "MasonryMacro2D::selectTrialFailureMode - element " << this->getTag() << " requires MasonryBendingMat or MasonryShearMat for bending and MasonryShearMat for shear" << endln;
        return -1;
    }

    const double shearRatio = std::abs(shearMaterial->getStrain()) / shearYieldStrain;
    const double bendingRatioI = std::abs(bendingMaterials[0]->getStrain()) / bendingYieldStrainI;
    const double bendingRatioJ = std::abs(bendingMaterials[1]->getStrain()) / bendingYieldStrainJ;
    const double bendingRatio = std::max(bendingRatioI, bendingRatioJ);
    if (shearRatio <= 1.0 && bendingRatio <= 1.0) {
        return 0;
    }
    trialFailureMode = shearRatio >= bendingRatio ? 1 : 2;
    return 0;
}

int MasonryMacro2D::solveInternalShearDeformationByNewton(double thetaI, double thetaJ, double axialForce, double initialShearDeformation)
{
    const double length = theCoordTransf->getInitialLength();
    double shearDeformation = initialShearDeformation;
    for (int iteration = 0; iteration < maximumIterations; ++iteration) {
        double residual = 0.0;
        double residualScale = 1.0;
        const int result = this->evaluateInternalShearResidual(thetaI, thetaJ, axialForce, shearDeformation, residual, residualScale);
        if (result != 0) {
            return result;
        }
        if (!useDisplacementConvergence && std::abs(residual) <= relativeTolerance * residualScale) {
            return 0;
        }

        const double tangentI = bendingMaterials[0]->getTangent();
        const double tangentJ = bendingMaterials[1]->getTangent();
        const double shearTangent = shearMaterial->getTangent();
        const double residualTangent = tangentI / length + tangentJ / length + shearTangent * length - axialForce;
        const double tangentScale = std::max(1.0, std::abs(tangentI / length) + std::abs(tangentJ / length) + std::abs(shearTangent * length) + std::abs(axialForce));
        if (std::abs(residualTangent) <= 100.0 * std::numeric_limits<double>::epsilon() * tangentScale) {
            break;
        }

        const double correction = residual / residualTangent;
        if (!std::isfinite(correction)) {
            break;
        }
        if (useDisplacementConvergence && std::abs(correction) <= relativeTolerance * std::max(1.0, std::abs(shearDeformation))) {
            return 0;
        }
        shearDeformation -= correction;
    }
    return -1;
}

int MasonryMacro2D::solveInternalShearDeformationByBisection(double thetaI, double thetaJ, double axialForce, double initialShearDeformation)
{
    const double length = theCoordTransf->getInitialLength();

    double centerResidual = 0.0;
    double centerResidualScale = 1.0;
    int result = this->evaluateInternalShearResidual(thetaI, thetaJ, axialForce, initialShearDeformation, centerResidual, centerResidualScale);
    if (result != 0) {
        return result;
    }
    if (!useDisplacementConvergence && std::abs(centerResidual) <= relativeTolerance * centerResidualScale) {
        return 0;
    }

    // residualA: 区间 A 端残差。
    // residualB: 区间 B 端残差。
    const auto hasSignChange = [](double residualA, double residualB) {
        return (residualA <= 0.0 && residualB >= 0.0) || (residualA >= 0.0 && residualB <= 0.0);
    };

    double lowerDeformation = initialShearDeformation;
    double upperDeformation = initialShearDeformation;
    double lowerResidual = centerResidual;
    double upperResidual = centerResidual;
    bool bracketFound = false;
    double searchRadius = std::max(1.0, length) * 1.0e-6;
    for (int search = 0; search < 50; ++search) {
        const double leftDeformation = initialShearDeformation - searchRadius;
        const double rightDeformation = initialShearDeformation + searchRadius;
        double leftResidual = 0.0;
        double rightResidual = 0.0;
        double trialResidualScale = 1.0;
        result = this->evaluateInternalShearResidual(thetaI, thetaJ, axialForce, leftDeformation, leftResidual, trialResidualScale);
        result += this->evaluateInternalShearResidual(thetaI, thetaJ, axialForce, rightDeformation, rightResidual, trialResidualScale);
        if (result != 0) {
            return result;
        }

        const bool leftBracket = hasSignChange(leftResidual, centerResidual);
        const bool rightBracket = hasSignChange(centerResidual, rightResidual);
        if (leftBracket || rightBracket) {
            const double leftEstimatedDistance = leftBracket ? searchRadius * std::abs(centerResidual) / (std::abs(leftResidual) + std::abs(centerResidual)) : std::numeric_limits<double>::max();
            const double rightEstimatedDistance = rightBracket ? searchRadius * std::abs(centerResidual) / (std::abs(rightResidual) + std::abs(centerResidual)) : std::numeric_limits<double>::max();
            if (leftEstimatedDistance <= rightEstimatedDistance) {
                lowerDeformation = leftDeformation;
                lowerResidual = leftResidual;
                upperDeformation = initialShearDeformation;
                upperResidual = centerResidual;
            } else {
                lowerDeformation = initialShearDeformation;
                lowerResidual = centerResidual;
                upperDeformation = rightDeformation;
                upperResidual = rightResidual;
            }
            bracketFound = true;
            break;
        }
        searchRadius *= 2.0;
    }

    // 在异号区间内二分，收敛点最后一次材料试算状态即为单元采用的横向状态。
    if (bracketFound) {
        for (int iteration = 0; iteration < maximumIterations; ++iteration) {
            const double middleDeformation = 0.5 * (lowerDeformation + upperDeformation);
            double middleResidual = 0.0;
            double middleResidualScale = 1.0;
            result = this->evaluateInternalShearResidual(thetaI, thetaJ, axialForce, middleDeformation, middleResidual, middleResidualScale);
            if (result != 0) {
                return result;
            }
            if (!useDisplacementConvergence && std::abs(middleResidual) <= relativeTolerance * middleResidualScale) {
                return 0;
            }

            if (hasSignChange(lowerResidual, middleResidual)) {
                upperDeformation = middleDeformation;
                upperResidual = middleResidual;
            } else {
                lowerDeformation = middleDeformation;
                lowerResidual = middleResidual;
            }
            if (useDisplacementConvergence && std::abs(upperDeformation - lowerDeformation) <= relativeTolerance * std::max(1.0, std::abs(middleDeformation))) {
                const double convergedDeformation = std::abs(lowerResidual) <= std::abs(upperResidual) ? lowerDeformation : upperDeformation;
                return this->evaluateInternalShearResidual(thetaI, thetaJ, axialForce, convergedDeformation, middleResidual, middleResidualScale);
            }
        }
    }
    return -1;
}

const Matrix &MasonryMacro2D::getTangentStiff(void)
//  用当前材料切线组装基本切线刚度，并通过 coordinateTransformation 转换到 6x6 全局刚度。
{
    // Mj = Mj(phiJ) = Mj(thetaJ + Delta / L)
    // Mi = Mi(phiI) = Mi(thetaI + Delta / L)
    // 
    // ∂Mj/∂thetaI = Mj' * ∂Delta/∂thetaI / L
    // ∂Mj/∂thetaJ = Mj' * (1 + ∂Delta/∂thetaJ / L) 
    // ∂Mi/∂thetaI = Mi' * (1 + ∂Delta/∂thetaI / L)
    // ∂Mi/∂thetaJ = Mi' * ∂Delta/∂thetaJ / L
    // 
    // 其中 ∂Delta/∂thetaI 和 ∂Delta/∂thetaJ 由平衡方程 Mi(phiI) + Mj(phiJ) + V(Delta) * L =  N * Delta 对 Delta 求导得到：
    // ∂Delta/∂thetaI = -Mi' / (Mi'/L + Mj'/L + V'*L - N)
    // ∂Delta/∂thetaJ = -Mj' / (Mi'/L + Mj'/L + V'*L - N)
    // 
    // 代入上式得到基本切线刚度矩阵：
    // ∂Mj/∂thetaI = -Mj' * Mi' / (Mi' + Mj' + V'*L^2 - N*L)
    // ∂Mj/∂thetaJ = Mj' * (1 - Mj' / (Mi' + Mj' + V'*L^2 - N*L))
    // ∂Mi/∂thetaI = Mi' * (1 - Mi' / (Mi' + Mj' + V'*L^2 - N*L))
    // ∂Mi/∂thetaJ = -Mi' * Mj' / (Mi' + Mj' + V'*L^2 - N*L)

    Matrix basicStiffness(3, 3);
    basicStiffness.Zero();

    // 轴向
    basicStiffness(0, 0) = axialMaterial->getTangent();

    // 横向切线采用冻结当前轴力的近似，忽略轴向变形引起的弯矩增量。
    const double Mj_prime = bendingMaterials[1]->getTangent();
    const double Mi_prime = bendingMaterials[0]->getTangent();
    const double V_prime = shearMaterial->getTangent();
    const double L = theCoordTransf->getInitialLength();
    const double N = axialMaterial->getStress();
    const double condensedDenominator = Mi_prime + Mj_prime + V_prime * L * L - N * L;
    const double denominatorScale = std::max(1.0, std::abs(Mi_prime) + std::abs(Mj_prime) + std::abs(V_prime * L * L) + std::abs(N * L));
    if (std::abs(condensedDenominator) <= 100.0 * std::numeric_limits<double>::epsilon() * denominatorScale) {
        opserr << "MasonryMacro2D::getTangentStiff - element " << this->getTag() << " has a singular condensed tangent" << endln;
        tangentStiffness.Zero();
        return tangentStiffness;
    }

    // 刚度
    basicStiffness(1, 1) = Mi_prime * (1.0 - Mi_prime / condensedDenominator);
    basicStiffness(1, 2) = -Mi_prime * Mj_prime / condensedDenominator;
    basicStiffness(2, 1) = basicStiffness(1, 2);
    basicStiffness(2, 2) = Mj_prime * (1.0 - Mj_prime / condensedDenominator);

    Vector basicForce(3);
    basicForce(0) = axialMaterial->getStress();
    basicForce(1) = bendingMaterials[0]->getStress();
    basicForce(2) = bendingMaterials[1]->getStress();

    // 通过坐标变换转换为全局刚度矩阵
    tangentStiffness = theCoordTransf->getGlobalStiffMatrix(basicStiffness, basicForce);

    return tangentStiffness;
}

const Matrix &MasonryMacro2D::getInitialStiff(void)
{
    Matrix initialBasicStiffness(3, 3);
    initialBasicStiffness.Zero();

    const double axialTangent = axialMaterial->getInitialTangent();
    const double tangentI = bendingMaterials[0]->getInitialTangent();
    const double tangentJ = bendingMaterials[1]->getInitialTangent();
    const double shearTangent = shearMaterial->getInitialTangent();
    const double length = theCoordTransf->getInitialLength();
    const double condensedDenominator = tangentI + tangentJ + shearTangent * length * length;
    const double denominatorScale = std::max(1.0, std::abs(tangentI) + std::abs(tangentJ) + std::abs(shearTangent * length * length));
    if (std::abs(condensedDenominator) <= 100.0 * std::numeric_limits<double>::epsilon() * denominatorScale) {
        opserr << "MasonryMacro2D::getInitialStiff - element " << this->getTag() << " has a singular initial condensed tangent" << endln;
        initialStiffness.Zero();
        return initialStiffness;
    }

    // 初始状态取 N=0，并对内部剪切变形执行与当前切线相同的静力凝聚。
    initialBasicStiffness(0, 0) = axialTangent;
    initialBasicStiffness(1, 1) = tangentI * (1.0 - tangentI / condensedDenominator);
    initialBasicStiffness(1, 2) = -tangentI * tangentJ / condensedDenominator;
    initialBasicStiffness(2, 1) = initialBasicStiffness(1, 2);
    initialBasicStiffness(2, 2) = tangentJ * (1.0 - tangentJ / condensedDenominator);

    initialStiffness = theCoordTransf->getInitialGlobalStiffMatrix(initialBasicStiffness);
    return initialStiffness;
}

const Vector &MasonryMacro2D::getResistingForce(void)
// 从各材料读取当前内力，组装并转换为 6 维全局节点恢复力。
{
    // 轴力
    double axialForce = axialMaterial->getStress();
    // Mi
    double momentI = bendingMaterials[0]->getStress();
    // Mj
    double momentJ = bendingMaterials[1]->getStress();

    Vector basicForce(3);
    basicForce (0) = axialForce;
    basicForce (1) = momentI;
    basicForce (2) = momentJ;

    Vector elementLoad(3);
    elementLoad.Zero();
    
    resistingForce = theCoordTransf->getGlobalResistingForce(basicForce, elementLoad);

    return resistingForce;
}

int MasonryMacro2D::sendSelf(int commitTag, Channel &theChannel)
{
    if (theCoordTransf == 0 || bendingMaterials[0] == 0 || bendingMaterials[1] == 0 || shearMaterial == 0 || axialMaterial == 0) {
        opserr << "MasonryMacro2D::sendSelf - element is not fully initialized" << endln;
        return -1;
    }

    int dataTag = this->getDbTag();
    if (theCoordTransf->getDbTag() == 0) {
        theCoordTransf->setDbTag(theChannel.getDbTag());
    }

    ID idData(13);
    idData(0) = this->getTag();
    idData(1) = connectedExternalNodes(0);
    idData(2) = connectedExternalNodes(1);
    idData(3) = theCoordTransf->getClassTag();
    idData(4) = theCoordTransf->getDbTag();
    UniaxialMaterial *materials[4] = {bendingMaterials[0], bendingMaterials[1], shearMaterial, axialMaterial};
    for (int i = 0; i < 4; ++i) {
        idData(5 + 2 * i) = materials[i]->getClassTag();
        idData(6 + 2 * i) = ensureMaterialDbTag(materials[i], theChannel);
    }

    Vector vectorData(17);
    vectorData(0) = alphaM;
    vectorData(1) = betaK;
    vectorData(2) = betaK0;
    vectorData(3) = betaKc;
    vectorData(4) = width;
    vectorData(5) = thickness;
    vectorData(6) = compressiveStrength;
    vectorData(7) = cohesion;
    vectorData(8) = diagonalTensileStrength;
    vectorData(9) = stressBlockCoefficient;
    vectorData(10) = frictionCoefficient;
    vectorData(11) = maximumIterations;
    vectorData(12) = relativeTolerance;
    vectorData(13) = useDisplacementConvergence ? 1.0 : 0.0;
    vectorData(14) = committedContraflexureDistance;
    vectorData(15) = useExclusiveFailureMode ? 1.0 : 0.0;
    vectorData(16) = committedFailureMode;

    if (theChannel.sendID(dataTag, commitTag, idData) < 0 || theChannel.sendVector(dataTag, commitTag, vectorData) < 0 || theCoordTransf->sendSelf(commitTag, theChannel) < 0) {
        return -1;
    }
    for (int i = 0; i < 4; ++i) {
        if (materials[i]->sendSelf(commitTag, theChannel) < 0) {
            return -1;
        }
    }
    // TODO: 若核心模型新增单元级历史变量，应同步扩展这里和 recvSelf 的数据布局。
    return 0;
}

int MasonryMacro2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    ID idData(13);
    Vector vectorData(17);
    if (theChannel.recvID(dataTag, commitTag, idData) < 0 || theChannel.recvVector(dataTag, commitTag, vectorData) < 0) {
        return -1;
    }

    this->setTag(idData(0));
    connectedExternalNodes(0) = idData(1);
    connectedExternalNodes(1) = idData(2);
    alphaM = vectorData(0);
    betaK = vectorData(1);
    betaK0 = vectorData(2);
    betaKc = vectorData(3);
    width = vectorData(4);
    thickness = vectorData(5);
    compressiveStrength = vectorData(6);
    cohesion = vectorData(7);
    diagonalTensileStrength = vectorData(8);
    stressBlockCoefficient = vectorData(9);
    frictionCoefficient = vectorData(10);
    maximumIterations = static_cast<int>(vectorData(11));
    relativeTolerance = vectorData(12);
    useDisplacementConvergence = vectorData(13) != 0.0;
    committedContraflexureDistance = vectorData(14);
    useExclusiveFailureMode = vectorData(15) != 0.0;
    committedFailureMode = static_cast<int>(vectorData(16));
    trialFailureMode = committedFailureMode;
    if (theCoordTransf == 0 || theCoordTransf->getClassTag() != idData(3)) {
        delete theCoordTransf;
        theCoordTransf = theBroker.getNewCrdTransf(idData(3));
    }
    if (theCoordTransf == 0) {
        return -1;
    }
    theCoordTransf->setDbTag(idData(4));
    if (theCoordTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
        return -1;
    }

    UniaxialMaterial **materials[4] = {&bendingMaterials[0], &bendingMaterials[1], &shearMaterial, &axialMaterial};
    for (int i = 0; i < 4; ++i) {
        if (receiveMaterial(*materials[i], idData(5 + 2 * i), idData(6 + 2 * i), commitTag, theChannel, theBroker) < 0) {
            return -1;
        }
    }
    // TODO: 若核心模型新增单元级历史变量，应同步扩展这里和 sendSelf 的数据布局。
    return 0;
}

int MasonryMacro2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numModes)
{
    if (theNodes[0] == 0 || theNodes[1] == 0) {
        return -1;
    }
    Vector pointI(3);
    Vector pointJ(3);
    theNodes[0]->getDisplayCrds(pointI, fact, displayMode);
    theNodes[1]->getDisplayCrds(pointJ, fact, displayMode);
    return theViewer.drawLine(pointI, pointJ, 1.0, 1.0, this->getTag());
}

void MasonryMacro2D::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        // 输出模型定义所需的标签、几何参数和局部迭代选项。
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"MasonryMacro2D\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"bendingMaterial\": \"" << bendingMaterials[0]->getTag() << "\", ";
        s << "\"shearMaterial\": \"" << shearMaterial->getTag() << "\", ";
        s << "\"axialMaterial\": \"" << axialMaterial->getTag() << "\", ";
        s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\", ";
        s << "\"axialForceInteraction\": " << (width != 0.0 ? "true" : "false") << ", ";
        s << "\"width\": " << width << ", ";
        s << "\"thickness\": " << thickness << ", ";
        s << "\"compressiveStrength\": " << compressiveStrength << ", ";
        s << "\"cohesion\": " << cohesion << ", ";
        s << "\"diagonalTensileStrength\": " << diagonalTensileStrength << ", ";
        s << "\"stressBlockCoefficient\": " << stressBlockCoefficient << ", ";
        s << "\"frictionCoefficient\": " << frictionCoefficient << ", ";
        s << "\"maximumIterations\": " << maximumIterations << ", ";
        s << "\"relativeTolerance\": " << relativeTolerance << ", ";
        s << "\"localIterationConvergence\": \"" << (useDisplacementConvergence ? "displacement" : "force") << "\", ";
        s << "\"exclusiveFailureMode\": " << (useExclusiveFailureMode ? "true" : "false") << "}";
        return;
    }

    s << "MasonryMacro2D, element: " << this->getTag() << ", nodes: " << connectedExternalNodes << endln;
    if (width != 0.0) {
        s << "  axial interaction: width=" << width << " thickness=" << thickness << " fm=" << compressiveStrength << " cohesion=" << cohesion << " ft=" << diagonalTensileStrength << " committedH0=" << committedContraflexureDistance << endln;
    }
    s << "  local iteration: maximumIterations=" << maximumIterations << " relativeTolerance=" << relativeTolerance << " convergence=" << (useDisplacementConvergence ? "displacement" : "force") << endln;
    if (useExclusiveFailureMode) {
        s << "  exclusive failure mode: " << committedFailureMode << " (0=unlocked, 1=shear, 2=flexural)" << endln;
    }
}

Response *MasonryMacro2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc == 0) {
        return 0;
    }
    if (std::strcmp(argv[0], "force") == 0 || std::strcmp(argv[0], "globalForce") == 0) {
        return new ElementResponse(this, 1, resistingForce);
    }
    if (std::strcmp(argv[0], "stiffness") == 0) {
        return new ElementResponse(this, 2, tangentStiffness);
    }
    if (std::strcmp(argv[0], "material") == 0 && argc > 2) {
        if (std::strcmp(argv[1], "bendingI") == 0) return bendingMaterials[0]->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "bendingJ") == 0) return bendingMaterials[1]->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "shear") == 0) return shearMaterial->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "axial") == 0) return axialMaterial->setResponse(&argv[2], argc - 2, output);
    }
    return this->Element::setResponse(argv, argc, output);
}

int MasonryMacro2D::getResponse(int responseID, Information &information)
{
    if (responseID == 1) return information.setVector(resistingForce);
    if (responseID == 2) return information.setMatrix(tangentStiffness);
    return this->Element::getResponse(responseID, information);
}
