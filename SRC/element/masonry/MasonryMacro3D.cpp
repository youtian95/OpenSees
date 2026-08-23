// 模块功能：实现三维砌体宏单元的非核心框架，包括命令解析、节点连接、状态管理、响应和序列化。
// 使用流程：先定义六类单轴材料和三维 geomTransf，再调用 element masonryMacro 创建本单元。
// 输入输出：输入节点、材料和坐标变换；核心刚度、恢复力与材料变形映射留在 TODO 中由开发者实现。

#include "MasonryMacro3D.h"

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

#include <cstring>

namespace {

// material: 需要分配数据库标签的材料。
// theChannel: 提供数据库标签的通讯通道。
int ensureMaterialDbTag3D(UniaxialMaterial *material, Channel &theChannel)
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
int receiveMaterial3D(UniaxialMaterial *&material, int classTag, int dbTag, int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
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

} // namespace

// 解析三维 masonryMacro 命令并创建单元。
// 命令格式：element masonryMacro tag iNode jNode bendingYMat bendingZMat shearYMat shearZMat torsionMat axialMat transfTag
void *OPS_MasonryMacro3D(void)
{
    if (OPS_GetNumRemainingInputArgs() != 10) {
        opserr << "WARNING incorrect arguments: element masonryMacro tag iNode jNode bendingYMat bendingZMat shearYMat shearZMat torsionMat axialMat transfTag" << endln;
        return 0;
    }

    int integerData[10];
    int numberOfData = 10;
    if (OPS_GetIntInput(&numberOfData, integerData) < 0) {
        opserr << "WARNING masonryMacro3D failed to read integer arguments" << endln;
        return 0;
    }

    UniaxialMaterial *bendingYMaterial = OPS_getUniaxialMaterial(integerData[3]);
    UniaxialMaterial *bendingZMaterial = OPS_getUniaxialMaterial(integerData[4]);
    UniaxialMaterial *shearYMaterial = OPS_getUniaxialMaterial(integerData[5]);
    UniaxialMaterial *shearZMaterial = OPS_getUniaxialMaterial(integerData[6]);
    UniaxialMaterial *torsionMaterial = OPS_getUniaxialMaterial(integerData[7]);
    UniaxialMaterial *axialMaterial = OPS_getUniaxialMaterial(integerData[8]);
    CrdTransf *coordinateTransformation = OPS_getCrdTransf(integerData[9]);

    if (bendingYMaterial == 0 || bendingZMaterial == 0 || shearYMaterial == 0 || shearZMaterial == 0 || torsionMaterial == 0 || axialMaterial == 0 || coordinateTransformation == 0) {
        opserr << "WARNING masonryMacro3D could not find a material or coordinate transformation for element " << integerData[0] << endln;
        return 0;
    }

    return new MasonryMacro3D(integerData[0], integerData[1], integerData[2], *bendingYMaterial, *bendingZMaterial, *shearYMaterial, *shearZMaterial, *torsionMaterial, *axialMaterial, *coordinateTransformation);
}

MasonryMacro3D::MasonryMacro3D(int tag, int nodeI, int nodeJ, UniaxialMaterial &bendingYMaterial, UniaxialMaterial &bendingZMaterial, UniaxialMaterial &shearYMaterialInput, UniaxialMaterial &shearZMaterialInput, UniaxialMaterial &torsionMaterialInput, UniaxialMaterial &axialMaterialInput, CrdTransf &coordinateTransformation)
    : Element(tag, ELE_TAG_MasonryMacro3D), connectedExternalNodes(2), theCoordTransf(coordinateTransformation.getCopy3d()), shearYMaterial(shearYMaterialInput.getCopy()), shearZMaterial(shearZMaterialInput.getCopy()), torsionMaterial(torsionMaterialInput.getCopy()), axialMaterial(axialMaterialInput.getCopy()), tangentStiffness(12, 12), initialStiffness(12, 12), resistingForce(12)
{
    connectedExternalNodes(0) = nodeI;
    connectedExternalNodes(1) = nodeJ;
    theNodes[0] = 0;
    theNodes[1] = 0;
    bendingYMaterials[0] = bendingYMaterial.getCopy();
    bendingYMaterials[1] = bendingYMaterial.getCopy();
    bendingZMaterials[0] = bendingZMaterial.getCopy();
    bendingZMaterials[1] = bendingZMaterial.getCopy();
}

MasonryMacro3D::MasonryMacro3D()
    : Element(0, ELE_TAG_MasonryMacro3D), connectedExternalNodes(2), theCoordTransf(0), shearYMaterial(0), shearZMaterial(0), torsionMaterial(0), axialMaterial(0), tangentStiffness(12, 12), initialStiffness(12, 12), resistingForce(12)
{
    theNodes[0] = 0;
    theNodes[1] = 0;
    bendingYMaterials[0] = 0;
    bendingYMaterials[1] = 0;
    bendingZMaterials[0] = 0;
    bendingZMaterials[1] = 0;
}

MasonryMacro3D::~MasonryMacro3D()
{
    delete theCoordTransf;
    delete bendingYMaterials[0];
    delete bendingYMaterials[1];
    delete bendingZMaterials[0];
    delete bendingZMaterials[1];
    delete shearYMaterial;
    delete shearZMaterial;
    delete torsionMaterial;
    delete axialMaterial;
}

int MasonryMacro3D::getNumExternalNodes(void) const
{
    return 2;
}

const ID &MasonryMacro3D::getExternalNodes(void)
{
    return connectedExternalNodes;
}

Node **MasonryMacro3D::getNodePtrs(void)
{
    return theNodes;
}

int MasonryMacro3D::getNumDOF(void)
{
    return 12;
}

void MasonryMacro3D::setDomain(Domain *theDomain)
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
        opserr << "MasonryMacro3D::setDomain - element " << this->getTag() << " could not find both nodes" << endln;
        return;
    }
    if (theNodes[0]->getNumberDOF() != 6 || theNodes[1]->getNumberDOF() != 6) {
        opserr << "MasonryMacro3D::setDomain - element " << this->getTag() << " requires six DOFs at each node" << endln;
        return;
    }
    if (theCoordTransf == 0 || theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
        opserr << "MasonryMacro3D::setDomain - element " << this->getTag() << " failed to initialize coordinate transformation" << endln;
        return;
    }

    if (theCoordTransf->getInitialLength() <= 0.0) {
        opserr << "MasonryMacro3D::setDomain - element " << this->getTag() << " has zero length" << endln;
        return;
    }
    this->DomainComponent::setDomain(theDomain);
}

int MasonryMacro3D::commitState(void)
{
    int result = this->Element::commitState();
    result += theCoordTransf->commitState();
    UniaxialMaterial *materials[8] = {bendingYMaterials[0], bendingYMaterials[1], bendingZMaterials[0], bendingZMaterials[1], shearYMaterial, shearZMaterial, torsionMaterial, axialMaterial};
    for (int i = 0; i < 8; ++i) result += materials[i]->commitState();
    return result;
}

int MasonryMacro3D::revertToLastCommit(void)
{
    int result = theCoordTransf->revertToLastCommit();
    UniaxialMaterial *materials[8] = {bendingYMaterials[0], bendingYMaterials[1], bendingZMaterials[0], bendingZMaterials[1], shearYMaterial, shearZMaterial, torsionMaterial, axialMaterial};
    for (int i = 0; i < 8; ++i) result += materials[i]->revertToLastCommit();
    return result;
}

int MasonryMacro3D::revertToStart(void)
{
    int result = theCoordTransf->revertToStart();
    UniaxialMaterial *materials[8] = {bendingYMaterials[0], bendingYMaterials[1], bendingZMaterials[0], bendingZMaterials[1], shearYMaterial, shearZMaterial, torsionMaterial, axialMaterial};
    for (int i = 0; i < 8; ++i) result += materials[i]->revertToStart();
    resistingForce.Zero();
    return result;
}

int MasonryMacro3D::update(void)
{
    int result = theCoordTransf->update();
    if (result != 0) {
        return result;
    }

    // TODO: 从坐标变换读取三维基本位移，计算轴向、双向剪切、扭转和两方向端部弯曲变形。
    // TODO: 将计算出的变形直接传给八个材料副本，明确两个弯曲轴和两个剪切方向的符号约定。
    return 0;
}

const Matrix &MasonryMacro3D::getTangentStiff(void)
{
    // TODO: 用当前材料切线组装三维基本切线刚度，并通过 coordinateTransformation 转换到 12x12 全局刚度。
    return tangentStiffness;
}

const Matrix &MasonryMacro3D::getInitialStiff(void)
{
    // TODO: 用各材料初始切线组装并缓存 12x12 初始全局刚度。
    return initialStiffness;
}

const Vector &MasonryMacro3D::getResistingForce(void)
{
    // TODO: 从各材料读取当前内力，组装并转换为 12 维全局节点恢复力。
    return resistingForce;
}

int MasonryMacro3D::sendSelf(int commitTag, Channel &theChannel)
{
    UniaxialMaterial *materials[8] = {bendingYMaterials[0], bendingYMaterials[1], bendingZMaterials[0], bendingZMaterials[1], shearYMaterial, shearZMaterial, torsionMaterial, axialMaterial};
    if (theCoordTransf == 0) {
        return -1;
    }
    for (int i = 0; i < 8; ++i) {
        if (materials[i] == 0) return -1;
    }

    int dataTag = this->getDbTag();
    if (theCoordTransf->getDbTag() == 0) {
        theCoordTransf->setDbTag(theChannel.getDbTag());
    }

    ID idData(21);
    idData(0) = this->getTag();
    idData(1) = connectedExternalNodes(0);
    idData(2) = connectedExternalNodes(1);
    idData(3) = theCoordTransf->getClassTag();
    idData(4) = theCoordTransf->getDbTag();
    for (int i = 0; i < 8; ++i) {
        idData(5 + 2 * i) = materials[i]->getClassTag();
        idData(6 + 2 * i) = ensureMaterialDbTag3D(materials[i], theChannel);
    }

    Vector vectorData(4);
    vectorData(0) = alphaM;
    vectorData(1) = betaK;
    vectorData(2) = betaK0;
    vectorData(3) = betaKc;

    if (theChannel.sendID(dataTag, commitTag, idData) < 0 || theChannel.sendVector(dataTag, commitTag, vectorData) < 0 || theCoordTransf->sendSelf(commitTag, theChannel) < 0) {
        return -1;
    }
    for (int i = 0; i < 8; ++i) {
        if (materials[i]->sendSelf(commitTag, theChannel) < 0) return -1;
    }
    // TODO: 若核心模型新增单元级历史变量，应同步扩展这里和 recvSelf 的数据布局。
    return 0;
}

int MasonryMacro3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    ID idData(21);
    Vector vectorData(4);
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
    if (theCoordTransf == 0 || theCoordTransf->getClassTag() != idData(3)) {
        delete theCoordTransf;
        theCoordTransf = theBroker.getNewCrdTransf(idData(3));
    }
    if (theCoordTransf == 0) return -1;
    theCoordTransf->setDbTag(idData(4));
    if (theCoordTransf->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;

    UniaxialMaterial **materials[8] = {&bendingYMaterials[0], &bendingYMaterials[1], &bendingZMaterials[0], &bendingZMaterials[1], &shearYMaterial, &shearZMaterial, &torsionMaterial, &axialMaterial};
    for (int i = 0; i < 8; ++i) {
        if (receiveMaterial3D(*materials[i], idData(5 + 2 * i), idData(6 + 2 * i), commitTag, theChannel, theBroker) < 0) return -1;
    }
    // TODO: 若核心模型新增单元级历史变量，应同步扩展这里和 sendSelf 的数据布局。
    return 0;
}

int MasonryMacro3D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numModes)
{
    if (theNodes[0] == 0 || theNodes[1] == 0) return -1;
    Vector pointI(3);
    Vector pointJ(3);
    theNodes[0]->getDisplayCrds(pointI, fact, displayMode);
    theNodes[1]->getDisplayCrds(pointJ, fact, displayMode);
    return theViewer.drawLine(pointI, pointJ, 1.0, 1.0, this->getTag());
}

void MasonryMacro3D::Print(OPS_Stream &s, int flag)
{
    s << "MasonryMacro3D, element: " << this->getTag() << ", nodes: " << connectedExternalNodes << endln;
}

Response *MasonryMacro3D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc == 0) return 0;
    if (std::strcmp(argv[0], "force") == 0 || std::strcmp(argv[0], "globalForce") == 0) return new ElementResponse(this, 1, resistingForce);
    if (std::strcmp(argv[0], "stiffness") == 0) return new ElementResponse(this, 2, tangentStiffness);
    if (std::strcmp(argv[0], "material") == 0 && argc > 2) {
        if (std::strcmp(argv[1], "bendingYI") == 0) return bendingYMaterials[0]->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "bendingYJ") == 0) return bendingYMaterials[1]->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "bendingZI") == 0) return bendingZMaterials[0]->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "bendingZJ") == 0) return bendingZMaterials[1]->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "shearY") == 0) return shearYMaterial->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "shearZ") == 0) return shearZMaterial->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "torsion") == 0) return torsionMaterial->setResponse(&argv[2], argc - 2, output);
        if (std::strcmp(argv[1], "axial") == 0) return axialMaterial->setResponse(&argv[2], argc - 2, output);
    }
    return this->Element::setResponse(argv, argc, output);
}

int MasonryMacro3D::getResponse(int responseID, Information &information)
{
    if (responseID == 1) return information.setVector(resistingForce);
    if (responseID == 2) return information.setMatrix(tangentStiffness);
    return this->Element::getResponse(responseID, information);
}
