/*
 * OpenSees uniaxial interface for Rinaldin et al. masonry springs.
 *
 * Usage: create RinaldinMasonryShear or RinaldinMasonryFlexural and
 * assign it to a translational or rotational zeroLength degree of freedom.
 * Inputs: 9 shear or 11 flexural parameters exported by So.ph.i.
 * Outputs: force/moment and tangent stiffness at the supplied deformation.
 */

#ifndef RinaldinMasonryMaterial_h
#define RinaldinMasonryMaterial_h

#include <UniaxialMaterial.h>

class Channel;
class FEM_ObjectBroker;
class OPS_Stream;

// Shared uniaxial implementation for the shear and flexural rules.
class RinaldinMasonryMaterial : public UniaxialMaterial
{
public:
    // Selects the backbone and unloading rules for each spring role.
    enum ModelType {
        Shear = 1,
        Flexural = 2
    };

    RinaldinMasonryMaterial(
        int tag,
        ModelType modelType,
        const double *parameters,
        int parameterCount);
    RinaldinMasonryMaterial();
    ~RinaldinMasonryMaterial() override = default;

    // Returns the type exposed by the corresponding Tcl command.
    const char *getClassType() const override;

    // Updates trial deformation, force, and tangent stiffness.
    int setTrialStrain(double strain, double strainRate = 0.0) override;
    double getStrain() override;
    double getStress() override;
    double getTangent() override;
    double getInitialTangent() override;

    // Commits, restores, or clears the material history.
    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;

    // Implements OpenSees copy, serialization, and output interfaces.
    UniaxialMaterial *getCopy() override;
    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(
        int commitTag,
        Channel &theChannel,
        FEM_ObjectBroker &theBroker) override;
    void Print(OPS_Stream &stream, int flag = 0) override;

private:
    static constexpr int MaxParameterCount = 11;

    // Stores extrema, energy, and piecewise path state for one step.
    struct State {
        double strain;
        double stress;
        double tangent;
        double maxPositiveStrain;
        double minNegativeStrain;
        double work;
        int direction;
        int pathStage;
        double pathStartStrain;
        double pathStartStress;
        double break1Strain;
        double break1Stress;
        double break2Strain;
        double break2Stress;
        double targetStrain;
        double targetStress;
    };

    ModelType modelType;
    int parameterCount;
    double parameters[MaxParameterCount];

    double elasticStiffness;
    double yieldForce;
    double firstPostYieldStiffness;
    double maximumForce;
    double secondPostYieldStiffness;
    double ultimateDeformation;
    double yieldDeformation;
    double peakDeformation;

    State committed;
    State trial;

    // Computes yield and peak deformation from the So.ph.i. inputs.
    void setDerivedParameters();
    bool validateParameters(bool printMessage) const;

    double envelopeStress(double strain) const;
    double envelopeTangent(double strain) const;
    double unloadingStiffness(double maximumAbsoluteStrain) const;
    double dissipatedEnergy(const State &state) const;

    // Creates piecewise unloading points at load reversal.
    void startReversalPath(int newDirection);
    void startShearReversalPath(int newDirection);
    void startFlexuralReversalPath(int newDirection);
    void evaluateCurrentPath(double strain);
    bool reached(double value, double boundary, int direction) const;
    static double lineValue(
        double strain,
        double strain1,
        double stress1,
        double strain2,
        double stress2,
        double fallbackStiffness,
        double &tangent);
};

// Parses shear Tcl arguments and returns an OpenSees material.
void *OPS_RinaldinMasonryShear();
// Parses flexural Tcl arguments and returns an OpenSees material.
void *OPS_RinaldinMasonryFlexural();

#endif
