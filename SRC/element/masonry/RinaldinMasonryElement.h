/*
 * Rinaldin masonry pier and spandrel macro-elements.
 *
 * The elements support two-dimensional ndm=2, ndf=3 and three-dimensional
 * ndm=3, ndf=6 models. Each element connects two external nodes and two
 * internal rotation nodes so the internal rotations remain global DOFs.
 * In 3D, only the wall-plane axial, shear, and rotation DOFs have stiffness.
 * Each element contains one axial spring, one shear spring, and two
 * independent flexural springs. The in-plane shear deformation is coupled
 * with the end rotations, and pier capacities are updated from the current
 * compressive axial force at every trial state.
 *
 * Tcl usage:
 *   element RinaldinMasonryPierElement tag iNode jNode axialMat shearMat flexMat
 *       width thickness height h0 fcm cohesion ft kd mu
 *   element RinaldinMasonrySpandrelElement tag iNode jNode axialMat shearMat flexMat
 *       depth thickness length cohesion ft interlock course alpha adjacentStress
 *       horizontalStress beta
 *
 * The material tags must refer to RinaldinMasonryMaterial-compatible uniaxial
 * materials. Units are inherited from the supplied materials and geometry.
 */

#ifndef RinaldinMasonryElement_h
#define RinaldinMasonryElement_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class FEM_ObjectBroker;
class Information;
class Node;
class OPS_Stream;
class Response;
class Renderer;
class UniaxialMaterial;

class RinaldinMasonryElement : public Element
{
public:
    RinaldinMasonryElement();
    RinaldinMasonryElement(int tag, int classTag, bool pier,
        int iExternal, int iInternal, int jInternal, int jExternal,
        UniaxialMaterial &axial, UniaxialMaterial &shear,
        UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
        UniaxialMaterial &rigidRotation,
        const double *geometry, int geometryCount);
    ~RinaldinMasonryElement() override;

    const char *getClassType() const override;
    int getNumExternalNodes() const override;
    const ID &getExternalNodes() override;
    Node **getNodePtrs() override;
    int getNumDOF() override;
    void setDomain(Domain *theDomain) override;

    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    int update() override;

    const Matrix &getTangentStiff() override;
    const Matrix &getInitialStiff() override;
    const Vector &getResistingForce() override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel,
        FEM_ObjectBroker &theBroker) override;
    void Print(OPS_Stream &s, int flag = 0) override;
    int displaySelf(Renderer &theViewer, int displayMode, float fact,
        const char **modes = 0, int numModes = 0) override;

    Response *setResponse(const char **argv, int argc,
        OPS_Stream &output) override;
    int getResponse(int responseID, Information &info) override;

protected:
    void initialize(int tag, int classTag, bool pier,
        int iExternal, int iInternal, int jInternal, int jExternal,
        UniaxialMaterial &axial, UniaxialMaterial &shear,
        UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
        UniaxialMaterial &rigidRotation,
        const double *geometry, int geometryCount);
    int formTrialState();
    void updateInteraction(double axialForce);
    void formTangent(bool initial);
    void transformLocalToGlobal(const Vector &local, Vector &global) const;
    void transformMatrixToGlobal(const Matrix &local, Matrix &global) const;

    bool isPier;
    bool isThreeDimensional;
    int numberDOF;
    ID connectedExternalNodes;
    Node *theNodes[4];
    UniaxialMaterial *axialMaterial[3];
    UniaxialMaterial *shearMaterial;
    UniaxialMaterial *flexuralMaterialI;
    UniaxialMaterial *flexuralMaterialJ;
    UniaxialMaterial *rigidTranslationMaterial;
    UniaxialMaterial *rigidRotationMaterial;

    double geometry[16];
    int geometryCount;
    double length;
    double cosine;
    double sine;
    double localAxis[3];
    double localInPlane[3];
    double localNormal[3];
    double baseShearMaximum;
    double baseFlexuralMaximum;
    double shearYieldRatio;
    double shearK1Ratio;
    double shearK2Ratio;
    double flexuralYieldRatio;
    double flexuralK1Ratio;
    double flexuralK2Ratio;

    Vector deformation;
    Vector localForce;
    Vector resistingForce;
    Matrix tangent;
    Matrix initialTangent;
};

class RinaldinMasonryPierElement : public RinaldinMasonryElement
{
public:
    RinaldinMasonryPierElement();
    RinaldinMasonryPierElement(int tag,
        int iExternal, int iInternal, int jInternal, int jExternal,
        UniaxialMaterial &axial, UniaxialMaterial &shear,
        UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
        UniaxialMaterial &rigidRotation,
        const double *geometry, int geometryCount = 9);
};

class RinaldinMasonrySpandrelElement : public RinaldinMasonryElement
{
public:
    RinaldinMasonrySpandrelElement();
    RinaldinMasonrySpandrelElement(int tag,
        int iExternal, int iInternal, int jInternal, int jExternal,
        UniaxialMaterial &axial, UniaxialMaterial &shear,
        UniaxialMaterial &flex, UniaxialMaterial &rigidTranslation,
        UniaxialMaterial &rigidRotation,
        const double *geometry, int geometryCount = 11);
};

void *OPS_RinaldinMasonryPierElement();
void *OPS_RinaldinMasonrySpandrelElement();

#endif
