/*
 * Implementation of the Rinaldin masonry pier and spandrel macro-elements.
 *
 * Each element connects two external nodes and two internal rotation nodes.
 * The internal rotations are explicit global DOFs, so the global solver
 * iterates on them instead of a hidden local Newton loop. The pier shear and
 * flexural capacities follow the trial axial force continuously inside
 * formTrialState; updating them only at commit made the residual jump
 * between steps and stalled Newton iterations under axial redistribution.
 */

#include "RinaldinMasonryElement.h"

#include "RinaldinMasonryMaterial.h"
#include <Channel.h>
#include <Domain.h>
#include <ElementResponse.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <Parameter.h>
#include <Renderer.h>
#include <UniaxialMaterial.h>
#include <classTags.h>
#include <elementAPI.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
constexpr double Tiny = 1.0e-12;

double positivePart(double value)
{
    return value > 0.0 ? value : 0.0;
}

double ratioOrZero(double numerator, double denominator)
{
    return std::abs(denominator) > Tiny ? numerator / denominator : 0.0;
}

void printPierUsage()
{
    opserr << "Want: element RinaldinMasonryPierElement tag iExternal "
           << "iInternal jInternal jExternal axialMat shearMat flexMat "
           << "rigidTransMat rigidRotMat width thickness height h0 fcm "
           << "cohesion ft kd mu [inPlaneX inPlaneY inPlaneZ]" << endln;
}

void printSpandrelUsage()
{
    opserr << "Want: element RinaldinMasonrySpandrelElement tag iExternal "
           << "iInternal jInternal jExternal axialMat shearMat flexMat "
           << "rigidTransMat rigidRotMat depth thickness length cohesion ft "
           << "interlock course alpha adjacentStress horizontalStress beta"
           << " [inPlaneX inPlaneY inPlaneZ]"
           << endln;
}

void updateMaterialEnvelope(UniaxialMaterial *material, double maximumForce,
    double yieldRatio, double k1Ratio, double k2Ratio)
{
    auto *rinaldin = dynamic_cast<RinaldinMasonryMaterial *>(material);
    if (rinaldin == nullptr) {
        return;
    }

    maximumForce = std::max(maximumForce, Tiny);
    const double yieldForce = std::max(yieldRatio * maximumForce, Tiny);
    const double k1 = k1Ratio * maximumForce;
    const double k2 = k2Ratio * maximumForce;
    rinaldin->setDynamicEnvelope(yieldForce, k1, maximumForce, k2);
}
}

void *OPS_RinaldinMasonryPierElement()
{
    const int numberArgs = OPS_GetNumRemainingInputArgs();
    if (numberArgs != 19 && numberArgs != 22) {
        printPierUsage();
        return nullptr;
    }

    int iData[10];
    int count = 10;
    if (OPS_GetIntInput(&count, iData) != 0) {
        opserr << "WARNING invalid integer data for Rinaldin pier element"
               << endln;
        return nullptr;
    }
    double geometry[12] = {};
    count = numberArgs - 10;
    if (OPS_GetDoubleInput(&count, geometry) != 0) {
        opserr << "WARNING invalid geometry data for Rinaldin pier element"
               << endln;
        return nullptr;
    }

    UniaxialMaterial *axial = OPS_getUniaxialMaterial(iData[5]);
    UniaxialMaterial *shear = OPS_getUniaxialMaterial(iData[6]);
    UniaxialMaterial *flex = OPS_getUniaxialMaterial(iData[7]);
    UniaxialMaterial *rigidTranslation = OPS_getUniaxialMaterial(iData[8]);
    UniaxialMaterial *rigidRotation = OPS_getUniaxialMaterial(iData[9]);
    if (axial == nullptr || shear == nullptr || flex == nullptr
        || rigidTranslation == nullptr || rigidRotation == nullptr) {
        opserr << "WARNING Rinaldin pier element material was not found"
               << endln;
        return nullptr;
    }
    return new RinaldinMasonryPierElement(
        iData[0], iData[1], iData[2], iData[3], iData[4],
        *axial, *shear, *flex, *rigidTranslation, *rigidRotation,
        geometry, count);
}

void *OPS_RinaldinMasonrySpandrelElement()
{
    const int numberArgs = OPS_GetNumRemainingInputArgs();
    if (numberArgs != 21 && numberArgs != 24) {
        printSpandrelUsage();
        return nullptr;
    }

    int iData[10];
    int count = 10;
    if (OPS_GetIntInput(&count, iData) != 0) {
        opserr << "WARNING invalid integer data for Rinaldin spandrel element"
               << endln;
        return nullptr;
    }
    double geometry[14] = {};
    count = numberArgs - 10;
    if (OPS_GetDoubleInput(&count, geometry) != 0) {
        opserr << "WARNING invalid geometry data for Rinaldin spandrel element"
               << endln;
        return nullptr;
    }

    UniaxialMaterial *axial = OPS_getUniaxialMaterial(iData[5]);
    UniaxialMaterial *shear = OPS_getUniaxialMaterial(iData[6]);
    UniaxialMaterial *flex = OPS_getUniaxialMaterial(iData[7]);
    UniaxialMaterial *rigidTranslation = OPS_getUniaxialMaterial(iData[8]);
    UniaxialMaterial *rigidRotation = OPS_getUniaxialMaterial(iData[9]);
    if (axial == nullptr || shear == nullptr || flex == nullptr
        || rigidTranslation == nullptr || rigidRotation == nullptr) {
        opserr << "WARNING Rinaldin spandrel element material was not found"
               << endln;
        return nullptr;
    }
    return new RinaldinMasonrySpandrelElement(
        iData[0], iData[1], iData[2], iData[3], iData[4],
        *axial, *shear, *flex, *rigidTranslation, *rigidRotation,
        geometry, count);
}

RinaldinMasonryElement::RinaldinMasonryElement()
    : Element(0, ELE_TAG_RinaldinMasonryPierElement),
      isPier(true), isThreeDimensional(false), numberDOF(6),
      connectedExternalNodes(4), theNodes{nullptr, nullptr, nullptr, nullptr},
      axialMaterial{nullptr, nullptr, nullptr}, shearMaterial(nullptr),
      flexuralMaterialI(nullptr), flexuralMaterialJ(nullptr),
      rigidTranslationMaterial(nullptr), rigidRotationMaterial(nullptr),
      geometryCount(0),
      length(0.0), cosine(1.0), sine(0.0), localAxis{1.0, 0.0, 0.0},
      localInPlane{0.0, 1.0, 0.0}, localNormal{0.0, 0.0, 1.0},
      baseShearMaximum(0.0),
      baseFlexuralMaximum(0.0), shearYieldRatio(0.0), shearK1Ratio(0.0),
      shearK2Ratio(0.0), flexuralYieldRatio(0.0), flexuralK1Ratio(0.0),
      flexuralK2Ratio(0.0), deformation(4), localForce(6), resistingForce(6),
      tangent(6, 6), initialTangent(6, 6)
{
    std::fill(geometry, geometry + 16, 0.0);
}

void RinaldinMasonryElement::initialize(int tag, int classTag, bool pier,
    int iExternal, int iInternal, int jInternal, int jExternal,
    UniaxialMaterial &axial, UniaxialMaterial &shear,
    UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
    UniaxialMaterial &rigidRotation,
    const double *values, int count)
{
    setTag(tag);
    isPier = pier;
    connectedExternalNodes(0) = iExternal;
    connectedExternalNodes(1) = iInternal;
    connectedExternalNodes(2) = jInternal;
    connectedExternalNodes(3) = jExternal;
    geometryCount = count;
    std::fill(geometry, geometry + 16, 0.0);
    std::copy(values, values + count, geometry);

    axialMaterial[0] = axial.getCopy();
    axialMaterial[1] = axial.getCopy();
    axialMaterial[2] = axial.getCopy();
    shearMaterial = shear.getCopy();
    flexuralMaterialI = flex.getCopy();
    flexuralMaterialJ = flex.getCopy();
    rigidTranslationMaterial = rigidTranslation.getCopy();
    rigidRotationMaterial = rigidRotation.getCopy();
    if (axialMaterial[0] == nullptr || axialMaterial[1] == nullptr
        || axialMaterial[2] == nullptr || shearMaterial == nullptr
        || flexuralMaterialI == nullptr || flexuralMaterialJ == nullptr
        || rigidTranslationMaterial == nullptr
        || rigidRotationMaterial == nullptr) {
        opserr << "FATAL Rinaldin masonry element failed to copy material"
               << endln;
        return;
    }

    baseShearMaximum = shearMaterial->getStress();
    baseFlexuralMaximum = flexuralMaterialI->getStress();
    auto *shearRinaldin = dynamic_cast<RinaldinMasonryMaterial *>(shearMaterial);
    auto *flexRinaldin = dynamic_cast<RinaldinMasonryMaterial *>(flexuralMaterialI);
    if (shearRinaldin != nullptr) {
        baseShearMaximum = shearRinaldin->getMaximumForce();
        shearYieldRatio = ratioOrZero(shearRinaldin->getYieldForce(), baseShearMaximum);
        shearK1Ratio = ratioOrZero(shearRinaldin->getFirstPostYieldStiffness(), baseShearMaximum);
        shearK2Ratio = ratioOrZero(shearRinaldin->getSecondPostYieldStiffness(), baseShearMaximum);
    }
    if (flexRinaldin != nullptr) {
        baseFlexuralMaximum = flexRinaldin->getMaximumForce();
        flexuralYieldRatio = ratioOrZero(flexRinaldin->getYieldForce(), baseFlexuralMaximum);
        flexuralK1Ratio = ratioOrZero(flexRinaldin->getFirstPostYieldStiffness(), baseFlexuralMaximum);
        flexuralK2Ratio = ratioOrZero(flexRinaldin->getSecondPostYieldStiffness(), baseFlexuralMaximum);
    }
}

RinaldinMasonryElement::RinaldinMasonryElement(int tag, int classTag, bool pier,
    int iExternal, int iInternal, int jInternal, int jExternal,
    UniaxialMaterial &axial, UniaxialMaterial &shear,
    UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
    UniaxialMaterial &rigidRotation,
    const double *values, int count)
    : Element(tag, classTag), isPier(pier),
      isThreeDimensional((pier && count >= 12) || (!pier && count >= 14)),
      numberDOF(isThreeDimensional ? 24 : 12), connectedExternalNodes(4),
      theNodes{nullptr, nullptr, nullptr, nullptr},
      axialMaterial{nullptr, nullptr, nullptr}, shearMaterial(nullptr),
      flexuralMaterialI(nullptr), flexuralMaterialJ(nullptr),
      rigidTranslationMaterial(nullptr), rigidRotationMaterial(nullptr),
      geometryCount(0),
      length(0.0), cosine(1.0), sine(0.0), localAxis{1.0, 0.0, 0.0},
      localInPlane{0.0, 1.0, 0.0}, localNormal{0.0, 0.0, 1.0},
      baseShearMaximum(0.0),
      baseFlexuralMaximum(0.0), shearYieldRatio(0.0), shearK1Ratio(0.0),
      shearK2Ratio(0.0), flexuralYieldRatio(0.0), flexuralK1Ratio(0.0),
      flexuralK2Ratio(0.0), deformation(4), localForce(numberDOF),
      resistingForce(numberDOF), tangent(numberDOF, numberDOF),
      initialTangent(numberDOF, numberDOF)
{
    std::fill(geometry, geometry + 16, 0.0);
    initialize(tag, classTag, pier, iExternal, iInternal, jInternal, jExternal,
        axial, shear, flex, rigidTranslation, rigidRotation, values, count);
}

RinaldinMasonryElement::~RinaldinMasonryElement()
{
    delete axialMaterial[0];
    delete axialMaterial[1];
    delete axialMaterial[2];
    delete shearMaterial;
    delete flexuralMaterialI;
    delete flexuralMaterialJ;
    delete rigidTranslationMaterial;
    delete rigidRotationMaterial;
}

const char *RinaldinMasonryElement::getClassType() const
{
    return isPier ? "RinaldinMasonryPierElement" : "RinaldinMasonrySpandrelElement";
}

int RinaldinMasonryElement::getNumExternalNodes() const { return 4; }
const ID &RinaldinMasonryElement::getExternalNodes() { return connectedExternalNodes; }
Node **RinaldinMasonryElement::getNodePtrs() { return theNodes; }
int RinaldinMasonryElement::getNumDOF() { return numberDOF; }

void RinaldinMasonryElement::setDomain(Domain *theDomain)
{
    for (int i = 0; i < 4; ++i) {
        theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
        if (theNodes[i] == nullptr) {
            opserr << "WARNING Rinaldin masonry element node was not found"
                   << endln;
            return;
        }
    }
    const int requiredDOF = isThreeDimensional ? 6 : 3;
    for (int i = 0; i < 4; ++i) {
        if (theNodes[i]->getNumberDOF() != requiredDOF) {
            opserr << "WARNING Rinaldin masonry element DOF count does not match "
                   << (isThreeDimensional ? "3D (ndf=6)" : "2D (ndf=3)")
                   << endln;
            return;
        }
    }

    const Vector &xi = theNodes[0]->getCrds();
    const Vector &xj = theNodes[3]->getCrds();
    double dx = xj(0) - xi(0);
    double dy = xj(1) - xi(1);
    double dz = isThreeDimensional ? xj(2) - xi(2) : 0.0;
    length = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (length <= Tiny) {
        opserr << "WARNING Rinaldin masonry element has zero length" << endln;
        return;
    }

    if (!isThreeDimensional) {
        cosine = dx / length;
        sine = dy / length;
        localAxis[0] = cosine;
        localAxis[1] = sine;
        localAxis[2] = 0.0;
        localInPlane[0] = -sine;
        localInPlane[1] = cosine;
        localInPlane[2] = 0.0;
        localNormal[0] = 0.0;
        localNormal[1] = 0.0;
        localNormal[2] = 1.0;
    } else {
        localAxis[0] = dx / length;
        localAxis[1] = dy / length;
        localAxis[2] = dz / length;
        const int orientationStart = isPier ? 9 : 11;
        if (geometryCount < orientationStart + 3) {
            opserr << "WARNING 3D Rinaldin masonry element requires an in-plane "
                   << "reference vector" << endln;
            return;
        }
        double reference[3] = {
            geometry[orientationStart], geometry[orientationStart + 1],
            geometry[orientationStart + 2]
        };
        const double projection = reference[0] * localAxis[0]
            + reference[1] * localAxis[1] + reference[2] * localAxis[2];
        for (int i = 0; i < 3; ++i) {
            localInPlane[i] = reference[i] - projection * localAxis[i];
        }
        const double inPlaneLength = std::sqrt(
            localInPlane[0] * localInPlane[0]
            + localInPlane[1] * localInPlane[1]
            + localInPlane[2] * localInPlane[2]);
        if (inPlaneLength <= Tiny) {
            opserr << "WARNING 3D in-plane reference is parallel to element axis"
                   << endln;
            return;
        }
        for (int i = 0; i < 3; ++i) {
            localInPlane[i] /= inPlaneLength;
        }
        localNormal[0] = localAxis[1] * localInPlane[2]
            - localAxis[2] * localInPlane[1];
        localNormal[1] = localAxis[2] * localInPlane[0]
            - localAxis[0] * localInPlane[2];
        localNormal[2] = localAxis[0] * localInPlane[1]
            - localAxis[1] * localInPlane[0];
    }

    DomainComponent::setDomain(theDomain);
    formTrialState();
    formTangent(true);
}

void RinaldinMasonryElement::updateInteraction(double axialForce)
{
    if (isPier) {
        const double b = geometry[0];
        const double t = geometry[1];
        const double h = geometry[2];
        const double h0 = geometry[3];
        const double fcm = geometry[4];
        const double cohesion = geometry[5];
        const double ft = geometry[6];
        const double kd = geometry[7];
        const double mu = geometry[8];
        const double compressionLimit = kd * b * t * fcm;
        const double n = std::min(std::max(axialForce, 0.0), compressionLimit);
        const double moment = n * b / 2.0 * (1.0 - n / compressionLimit);
        const double sliding = n > Tiny
            ? (1.5 * b * t * cohesion + mu * n)
                / (1.0 + 3.0 * h0 * t * cohesion / n)
            : 0.0;
        const double tensileArea = ft * b * t;
        const double xi = std::min(std::max(h / b, 1.0), 1.5);
        const double diagonal = tensileArea / xi
            * std::sqrt(1.0 + n / tensileArea);
        updateMaterialEnvelope(shearMaterial, std::min(sliding, diagonal),
            shearYieldRatio, shearK1Ratio, shearK2Ratio);
        updateMaterialEnvelope(flexuralMaterialI, moment, flexuralYieldRatio,
            flexuralK1Ratio, flexuralK2Ratio);
        updateMaterialEnvelope(flexuralMaterialJ, moment, flexuralYieldRatio,
            flexuralK1Ratio, flexuralK2Ratio);
    } else {
        const double depth = geometry[0];
        const double thickness = geometry[1];
        const double cohesion = geometry[3];
        const double ft = geometry[4];
        const double interlock = geometry[5];
        const double course = geometry[6];
        const double alpha = geometry[7];
        const double adjacentStress = geometry[8];
        const double horizontalStress = geometry[9];
        const double beta = geometry[10];
        const double equivalentTension = interlock / course
            * (cohesion + alpha * adjacentStress);
        const double xi = depth / (4.0 * course);
        const double moment = 2.0 / 3.0 * equivalentTension * thickness
            * interlock * depth * xi;
        const double shear = ft * depth * thickness * beta
            * std::sqrt(1.0 + horizontalStress / ft);
        updateMaterialEnvelope(shearMaterial, shear, shearYieldRatio,
            shearK1Ratio, shearK2Ratio);
        updateMaterialEnvelope(flexuralMaterialI, moment, flexuralYieldRatio,
            flexuralK1Ratio, flexuralK2Ratio);
        updateMaterialEnvelope(flexuralMaterialJ, moment, flexuralYieldRatio,
            flexuralK1Ratio, flexuralK2Ratio);
    }
}

int RinaldinMasonryElement::formTrialState()
{
    if (theNodes[0] == nullptr || theNodes[1] == nullptr
        || theNodes[2] == nullptr || theNodes[3] == nullptr) {
        return -1;
    }

    double nodeLocal[4][6] = {};
    for (int node = 0; node < 4; ++node) {
        const Vector &disp = theNodes[node]->getTrialDisp();
        if (!isThreeDimensional) {
            nodeLocal[node][0] = cosine * disp(0) + sine * disp(1);
            nodeLocal[node][1] = -sine * disp(0) + cosine * disp(1);
            nodeLocal[node][2] = disp(2);
        } else {
            nodeLocal[node][0] = localAxis[0] * disp(0)
                + localAxis[1] * disp(1) + localAxis[2] * disp(2);
            nodeLocal[node][1] = localInPlane[0] * disp(0)
                + localInPlane[1] * disp(1) + localInPlane[2] * disp(2);
            nodeLocal[node][2] = localNormal[0] * disp(0)
                + localNormal[1] * disp(1) + localNormal[2] * disp(2);
            nodeLocal[node][3] = localAxis[0] * disp(3)
                + localAxis[1] * disp(4) + localAxis[2] * disp(5);
            nodeLocal[node][4] = localInPlane[0] * disp(3)
                + localInPlane[1] * disp(4) + localInPlane[2] * disp(5);
            nodeLocal[node][5] = localNormal[0] * disp(3)
                + localNormal[1] * disp(4) + localNormal[2] * disp(5);
        }
    }

    // Update the interaction envelope with the trial axial force so the
    // resisting force stays continuous across committed steps.
    axialMaterial[1]->setTrialStrain(nodeLocal[2][0] - nodeLocal[1][0]);
    updateInteraction(std::abs(axialMaterial[1]->getStress()));

    localForce.Zero();
    tangent.Zero();
    deformation.Zero();

    auto addSpring = [&](int a, int b, bool central,
        UniaxialMaterial *mat0, UniaxialMaterial *mat1,
        UniaxialMaterial *mat2, UniaxialMaterial *mat3,
        UniaxialMaterial *mat4, UniaxialMaterial *mat5) {
        const int componentCount = isThreeDimensional ? 6 : 3;
        const int nodeDof = isThreeDimensional ? 6 : 3;
        double deformationValues[6] = {};
        double materialForce[6] = {};
        double materialTangent[6] = {};
        UniaxialMaterial *materials[6] = {
            mat0, mat1, mat2, mat3, mat4, mat5
        };

        deformationValues[0] = nodeLocal[b][0] - nodeLocal[a][0];
        deformationValues[1] = nodeLocal[b][1] - nodeLocal[a][1];
        const int rotationIndex = isThreeDimensional ? 5 : 2;
        if (central) {
            deformationValues[1] -= 0.5 * length
                * (nodeLocal[a][rotationIndex] + nodeLocal[b][rotationIndex]);
        }
        if (isThreeDimensional) {
            deformationValues[2] = nodeLocal[b][2] - nodeLocal[a][2];
            deformationValues[3] = nodeLocal[b][3] - nodeLocal[a][3];
            deformationValues[4] = nodeLocal[b][4] - nodeLocal[a][4];
            deformationValues[5] = nodeLocal[b][5] - nodeLocal[a][5];
        } else {
            deformationValues[2] = nodeLocal[b][2] - nodeLocal[a][2];
        }

        for (int i = 0; i < componentCount; ++i) {
            materials[i]->setTrialStrain(deformationValues[i]);
            materialForce[i] = materials[i]->getStress();
            materialTangent[i] = materials[i]->getTangent();
        }

        double B[6][12] = {};
        const int pairA = 0;
        const int pairB = nodeDof;
        B[0][pairA] = -1.0;
        B[0][pairB] = 1.0;
        B[1][pairA + 1] = -1.0;
        B[1][pairB + 1] = 1.0;
        if (central) {
            B[1][pairA + rotationIndex] = -0.5 * length;
            B[1][pairB + rotationIndex] = -0.5 * length;
        }
        B[2][pairA + 2] = -1.0;
        B[2][pairB + 2] = 1.0;
        if (isThreeDimensional) {
            B[3][pairA + 3] = -1.0;
            B[3][pairB + 3] = 1.0;
            B[4][pairA + 4] = -1.0;
            B[4][pairB + 4] = 1.0;
            B[5][pairA + 5] = -1.0;
            B[5][pairB + 5] = 1.0;
        }

        const int globalA = a * nodeDof;
        const int globalB = b * nodeDof;
        const int pairDof = 2 * nodeDof;
        for (int col = 0; col < pairDof; ++col) {
            const int node = col < nodeDof ? a : b;
            const int component = col % nodeDof;
            const int global = node * nodeDof + component;
            for (int k = 0; k < componentCount; ++k) {
                localForce(global) += B[k][col] * materialForce[k];
            }
        }

        for (int colI = 0; colI < pairDof; ++colI) {
            const int nodeI = colI < nodeDof ? a : b;
            const int componentI = colI % nodeDof;
            const int globalI = nodeI * nodeDof + componentI;
            for (int colJ = 0; colJ < pairDof; ++colJ) {
                const int nodeJ = colJ < nodeDof ? a : b;
                const int componentJ = colJ % nodeDof;
                const int globalJ = nodeJ * nodeDof + componentJ;
                for (int k = 0; k < componentCount; ++k) {
                    tangent(globalI, globalJ) +=
                        B[k][colI] * materialTangent[k] * B[k][colJ];
                }
            }
        }
    };

    if (isThreeDimensional) {
        addSpring(0, 1, false,
            axialMaterial[0], rigidTranslationMaterial,
            rigidTranslationMaterial, rigidRotationMaterial,
            rigidRotationMaterial, flexuralMaterialI);
        addSpring(1, 2, true,
            axialMaterial[1], shearMaterial, rigidTranslationMaterial,
            rigidRotationMaterial, rigidRotationMaterial,
            rigidRotationMaterial);
        addSpring(2, 3, false,
            axialMaterial[2], rigidTranslationMaterial,
            rigidTranslationMaterial, rigidRotationMaterial,
            rigidRotationMaterial, flexuralMaterialJ);
    } else {
        addSpring(0, 1, false,
            axialMaterial[0], rigidTranslationMaterial, flexuralMaterialI,
            nullptr, nullptr, nullptr);
        addSpring(1, 2, true,
            axialMaterial[1], shearMaterial, rigidRotationMaterial,
            nullptr, nullptr, nullptr);
        addSpring(2, 3, false,
            axialMaterial[2], rigidTranslationMaterial, flexuralMaterialJ,
            nullptr, nullptr, nullptr);
    }

    Matrix globalTangent(numberDOF, numberDOF);
    transformMatrixToGlobal(tangent, globalTangent);
    tangent = globalTangent;

    deformation(0) = nodeLocal[3][0] - nodeLocal[0][0];
    const int rotationIndex = isThreeDimensional ? 5 : 2;
    deformation(1) = (nodeLocal[2][1] - nodeLocal[1][1])
        - 0.5 * length * (nodeLocal[1][rotationIndex]
            + nodeLocal[2][rotationIndex]);
    deformation(2) = nodeLocal[1][2] - nodeLocal[0][2];
    deformation(3) = nodeLocal[3][2] - nodeLocal[2][2];

    transformLocalToGlobal(localForce, resistingForce);
    return 0;
}

int RinaldinMasonryElement::update()
{
    return formTrialState();
}

void RinaldinMasonryElement::formTangent(bool initial)
{
    Matrix local(numberDOF, numberDOF);
    local.Zero();
    const int rotationIndex = isThreeDimensional ? 5 : 2;

    auto addSpringTangent = [&](int a, int b, bool central,
        UniaxialMaterial *mat0, UniaxialMaterial *mat1,
        UniaxialMaterial *mat2, UniaxialMaterial *mat3,
        UniaxialMaterial *mat4, UniaxialMaterial *mat5) {
        const int componentCount = isThreeDimensional ? 6 : 3;
        const int nodeDof = isThreeDimensional ? 6 : 3;
        double materialTangent[6] = {};
        UniaxialMaterial *materials[6] = {
            mat0, mat1, mat2, mat3, mat4, mat5
        };
        for (int i = 0; i < componentCount; ++i) {
            materialTangent[i] = initial
                ? materials[i]->getInitialTangent()
                : materials[i]->getTangent();
        }

        double B[6][12] = {};
        const int pairA = 0;
        const int pairB = nodeDof;
        B[0][pairA] = -1.0;
        B[0][pairB] = 1.0;
        B[1][pairA + 1] = -1.0;
        B[1][pairB + 1] = 1.0;
        if (central) {
            B[1][pairA + rotationIndex] = -0.5 * length;
            B[1][pairB + rotationIndex] = -0.5 * length;
        }
        B[2][pairA + 2] = -1.0;
        B[2][pairB + 2] = 1.0;
        if (isThreeDimensional) {
            B[3][pairA + 3] = -1.0;
            B[3][pairB + 3] = 1.0;
            B[4][pairA + 4] = -1.0;
            B[4][pairB + 4] = 1.0;
            B[5][pairA + 5] = -1.0;
            B[5][pairB + 5] = 1.0;
        }

        const int globalA = a * nodeDof;
        const int globalB = b * nodeDof;
        const int pairDof = 2 * nodeDof;
        for (int colI = 0; colI < pairDof; ++colI) {
            const int nodeI = colI < nodeDof ? a : b;
            const int componentI = colI % nodeDof;
            const int globalI = nodeI * nodeDof + componentI;
            for (int colJ = 0; colJ < pairDof; ++colJ) {
                const int nodeJ = colJ < nodeDof ? a : b;
                const int componentJ = colJ % nodeDof;
                const int globalJ = nodeJ * nodeDof + componentJ;
                for (int k = 0; k < componentCount; ++k) {
                    local(globalI, globalJ) +=
                        B[k][colI] * materialTangent[k] * B[k][colJ];
                }
            }
        }
    };

    if (isThreeDimensional) {
        addSpringTangent(0, 1, false,
            axialMaterial[0], rigidTranslationMaterial,
            rigidTranslationMaterial, rigidRotationMaterial,
            rigidRotationMaterial, flexuralMaterialI);
        addSpringTangent(1, 2, true,
            axialMaterial[1], shearMaterial, rigidTranslationMaterial,
            rigidRotationMaterial, rigidRotationMaterial,
            rigidRotationMaterial);
        addSpringTangent(2, 3, false,
            axialMaterial[2], rigidTranslationMaterial,
            rigidTranslationMaterial, rigidRotationMaterial,
            rigidRotationMaterial, flexuralMaterialJ);
    } else {
        addSpringTangent(0, 1, false,
            axialMaterial[0], rigidTranslationMaterial, flexuralMaterialI,
            nullptr, nullptr, nullptr);
        addSpringTangent(1, 2, true,
            axialMaterial[1], shearMaterial, rigidRotationMaterial,
            nullptr, nullptr, nullptr);
        addSpringTangent(2, 3, false,
            axialMaterial[2], rigidTranslationMaterial, flexuralMaterialJ,
            nullptr, nullptr, nullptr);
    }

    if (initial) {
        transformMatrixToGlobal(local, initialTangent);
    } else {
        transformMatrixToGlobal(local, tangent);
    }
}

void RinaldinMasonryElement::transformLocalToGlobal(const Vector &local, Vector &global) const
{
    global.Zero();
    for (int node = 0; node < 4; ++node) {
        const int base = node * (isThreeDimensional ? 6 : 3);
        if (!isThreeDimensional) {
            global(base) = cosine * local(base) - sine * local(base + 1);
            global(base + 1) = sine * local(base) + cosine * local(base + 1);
            global(base + 2) = local(base + 2);
        } else {
            const double fx = local(base);
            const double fy = local(base + 1);
            const double fz = local(base + 2);
            const double mx = local(base + 3);
            const double my = local(base + 4);
            const double mz = local(base + 5);
            global(base) = localAxis[0] * fx + localInPlane[0] * fy
                + localNormal[0] * fz;
            global(base + 1) = localAxis[1] * fx + localInPlane[1] * fy
                + localNormal[1] * fz;
            global(base + 2) = localAxis[2] * fx + localInPlane[2] * fy
                + localNormal[2] * fz;
            global(base + 3) = localAxis[0] * mx + localInPlane[0] * my
                + localNormal[0] * mz;
            global(base + 4) = localAxis[1] * mx + localInPlane[1] * my
                + localNormal[1] * mz;
            global(base + 5) = localAxis[2] * mx + localInPlane[2] * my
                + localNormal[2] * mz;
        }
    }
}

void RinaldinMasonryElement::transformMatrixToGlobal(const Matrix &local, Matrix &global) const
{
    const int nodeDof = isThreeDimensional ? 6 : 3;
    double T[24][24] = {};
    if (!isThreeDimensional) {
        for (int node = 0; node < 4; ++node) {
            const int base = 3 * node;
            T[base][base] = cosine;
            T[base][base + 1] = sine;
            T[base + 1][base] = -sine;
            T[base + 1][base + 1] = cosine;
            T[base + 2][base + 2] = 1.0;
        }
    } else {
        const double axes[3][3] = {
            {localAxis[0], localAxis[1], localAxis[2]},
            {localInPlane[0], localInPlane[1], localInPlane[2]},
            {localNormal[0], localNormal[1], localNormal[2]}
        };
        for (int node = 0; node < 4; ++node) {
            const int base = 6 * node;
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < 3; ++col) {
                    T[base + row][base + col] = axes[row][col];
                    T[base + 3 + row][base + 3 + col] = axes[row][col];
                }
            }
        }
    }

    global.Zero();
    for (int i = 0; i < numberDOF; ++i) {
        for (int j = 0; j < numberDOF; ++j) {
            for (int a = 0; a < numberDOF; ++a) {
                for (int b = 0; b < numberDOF; ++b) {
                    global(i, j) += T[a][i] * local(a, b) * T[b][j];
                }
            }
        }
    }
}

const Matrix &RinaldinMasonryElement::getTangentStiff()
{
    formTrialState();
    return tangent;
}

const Matrix &RinaldinMasonryElement::getInitialStiff() { return initialTangent; }

const Vector &RinaldinMasonryElement::getResistingForce()
{
    formTrialState();
    return resistingForce;
}

int RinaldinMasonryElement::commitState()
{
    int result = Element::commitState();
    for (int i = 0; i < 3; ++i) {
        result += axialMaterial[i]->commitState();
    }
    result += shearMaterial->commitState();
    result += flexuralMaterialI->commitState();
    result += flexuralMaterialJ->commitState();
    result += rigidTranslationMaterial->commitState();
    result += rigidRotationMaterial->commitState();
    return result;
}

int RinaldinMasonryElement::revertToLastCommit()
{
    int result = 0;
    for (int i = 0; i < 3; ++i) {
        result += axialMaterial[i]->revertToLastCommit();
    }
    result += shearMaterial->revertToLastCommit();
    result += flexuralMaterialI->revertToLastCommit();
    result += flexuralMaterialJ->revertToLastCommit();
    result += rigidTranslationMaterial->revertToLastCommit();
    result += rigidRotationMaterial->revertToLastCommit();
    return result;
}

int RinaldinMasonryElement::revertToStart()
{
    int result = 0;
    for (int i = 0; i < 3; ++i) {
        result += axialMaterial[i]->revertToStart();
    }
    result += shearMaterial->revertToStart();
    result += flexuralMaterialI->revertToStart();
    result += flexuralMaterialJ->revertToStart();
    result += rigidTranslationMaterial->revertToStart();
    result += rigidRotationMaterial->revertToStart();
    deformation.Zero();
    localForce.Zero();
    resistingForce.Zero();
    return result;
}

int RinaldinMasonryElement::sendSelf(int, Channel &) { return -1; }
int RinaldinMasonryElement::recvSelf(int, Channel &, FEM_ObjectBroker &) { return -1; }

void RinaldinMasonryElement::Print(OPS_Stream &s, int)
{
    s << getClassType() << " tag: " << getTag()
      << " nodes: ";
    for (int i = 0; i < 4; ++i) {
        s << connectedExternalNodes(i) << " ";
    }
    s << endln;
}

int RinaldinMasonryElement::displaySelf(Renderer &theViewer, int displayMode,
    float fact, const char **, int)
{
    if (theNodes[0] == nullptr || theNodes[3] == nullptr) {
        return 0;
    }
    Vector xI(3), xJ(3);
    const Vector &cI = theNodes[0]->getCrds();
    const Vector &cJ = theNodes[3]->getCrds();
    const Vector &uI = theNodes[0]->getDisp();
    const Vector &uJ = theNodes[3]->getDisp();
    xI(0) = cI(0) + (displayMode ? fact * uI(0) : 0.0);
    xI(1) = cI(1) + (displayMode ? fact * uI(1) : 0.0);
    xI(2) = isThreeDimensional
        ? cI(2) + (displayMode ? fact * uI(2) : 0.0) : 0.0;
    xJ(0) = cJ(0) + (displayMode ? fact * uJ(0) : 0.0);
    xJ(1) = cJ(1) + (displayMode ? fact * uJ(1) : 0.0);
    xJ(2) = isThreeDimensional
        ? cJ(2) + (displayMode ? fact * uJ(2) : 0.0) : 0.0;
    return theViewer.drawLine(xI, xJ, 1.0, 1.0);
}

Response *RinaldinMasonryElement::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    if (argc == 0) {
        return nullptr;
    }
    if (std::strcmp(argv[0], "localForce") == 0) {
        return new ElementResponse(this, 1, localForce);
    }
    if (std::strcmp(argv[0], "deformation") == 0) {
        return new ElementResponse(this, 2, deformation);
    }
    if (std::strcmp(argv[0], "stiffness") == 0) {
        return new ElementResponse(this, 3, getTangentStiff());
    }
    return Element::setResponse(argv, argc, output);
}

int RinaldinMasonryElement::getResponse(int responseID, Information &info)
{
    formTrialState();
    if (responseID == 1) {
        return info.setVector(localForce);
    }
    if (responseID == 2) {
        return info.setVector(deformation);
    }
    if (responseID == 3) {
        return info.setMatrix(tangent);
    }
    return Element::getResponse(responseID, info);
}

RinaldinMasonryPierElement::RinaldinMasonryPierElement()
    : RinaldinMasonryElement()
{
    isPier = true;
}

RinaldinMasonryPierElement::RinaldinMasonryPierElement(int tag,
    int iExternal, int iInternal, int jInternal, int jExternal,
    UniaxialMaterial &axial, UniaxialMaterial &shear,
    UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
    UniaxialMaterial &rigidRotation,
    const double *values, int geometryCount)
    : RinaldinMasonryElement(tag, ELE_TAG_RinaldinMasonryPierElement, true,
        iExternal, iInternal, jInternal, jExternal,
        axial, shear, flex, rigidTranslation, rigidRotation,
        values, geometryCount)
{
}

RinaldinMasonrySpandrelElement::RinaldinMasonrySpandrelElement()
    : RinaldinMasonryElement()
{
    isPier = false;
}

RinaldinMasonrySpandrelElement::RinaldinMasonrySpandrelElement(int tag,
    int iExternal, int iInternal, int jInternal, int jExternal,
    UniaxialMaterial &axial, UniaxialMaterial &shear,
    UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
    UniaxialMaterial &rigidRotation,
    const double *values, int geometryCount)
    : RinaldinMasonryElement(tag, ELE_TAG_RinaldinMasonrySpandrelElement, false,
        iExternal, iInternal, jInternal, jExternal,
        axial, shear, flex, rigidTranslation, rigidRotation,
        values, geometryCount)
{
}
