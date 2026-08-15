/*
 * This module implements the trilinear shear hysteresis and three-stage
 * flexural unloading rules described by Rinaldin et al. (2016) and So.ph.i.
 *
 * Usage: Tcl factories read the tag and parameters. The material receives
 * zeroLength deformation and returns force/moment and tangent stiffness.
 * Input and output units must remain consistent with the calling model.
 */

#include "RinaldinMasonryMaterial.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
constexpr double Tiny = 1.0e-12;

double signOf(double value)
{
    return value < 0.0 ? -1.0 : 1.0;
}

void printShearUsage()
{
    opserr << "Want: uniaxialMaterial RinaldinMasonryShear tag "
           << "Kel Fy K1pl Fmax K2pl CF alpha beta Uult" << endln;
}

void printFlexuralUsage()
{
    opserr << "Want: uniaxialMaterial RinaldinMasonryFlexural tag "
           << "Kel Fy K1pl Fmax CC CF alpha CD Uult K2pl" << endln;
}
} // namespace

void *OPS_RinaldinMasonryShear()
{
    if (OPS_GetNumRemainingInputArgs() != 10) {
        printShearUsage();
        return nullptr;
    }

    int tag = 0;
    int count = 1;
    if (OPS_GetIntInput(&count, &tag) != 0) {
        opserr << "WARNING invalid tag for RinaldinMasonryShear" << endln;
        return nullptr;
    }

    double values[9];
    count = 9;
    if (OPS_GetDoubleInput(&count, values) != 0) {
        opserr << "WARNING invalid parameters for RinaldinMasonryShear "
               << tag << endln;
        return nullptr;
    }

    auto *material = new RinaldinMasonryMaterial(
        tag, RinaldinMasonryMaterial::Shear, values, 9);
    if (material == nullptr || material->getInitialTangent() <= 0.0) {
        delete material;
        return nullptr;
    }
    return material;
}

void *OPS_RinaldinMasonryFlexural()
{
    if (OPS_GetNumRemainingInputArgs() != 11) {
        printFlexuralUsage();
        return nullptr;
    }

    int tag = 0;
    int count = 1;
    if (OPS_GetIntInput(&count, &tag) != 0) {
        opserr << "WARNING invalid tag for RinaldinMasonryFlexural" << endln;
        return nullptr;
    }

    double values[10];
    count = 10;
    if (OPS_GetDoubleInput(&count, values) != 0) {
        opserr << "WARNING invalid parameters for RinaldinMasonryFlexural "
               << tag << endln;
        return nullptr;
    }

    auto *material = new RinaldinMasonryMaterial(
        tag, RinaldinMasonryMaterial::Flexural, values, 10);
    if (material == nullptr || material->getInitialTangent() <= 0.0) {
        delete material;
        return nullptr;
    }
    return material;
}

RinaldinMasonryMaterial::RinaldinMasonryMaterial(
    int tag,
    ModelType type,
    const double *values,
    int count)
    : UniaxialMaterial(tag, MAT_TAG_RinaldinMasonryMaterial),
      modelType(type),
      parameterCount(count),
      elasticStiffness(0.0),
      yieldForce(0.0),
      firstPostYieldStiffness(0.0),
      maximumForce(0.0),
      secondPostYieldStiffness(0.0),
      ultimateDeformation(0.0),
      yieldDeformation(0.0),
      peakDeformation(0.0)
{
    std::fill(parameters, parameters + MaxParameterCount, 0.0);
    if (values != nullptr && count > 0 && count <= MaxParameterCount) {
        std::copy(values, values + count, parameters);
    }
    setDerivedParameters();
    if (!validateParameters(true)) {
        elasticStiffness = -1.0;
    }
    revertToStart();
}

RinaldinMasonryMaterial::RinaldinMasonryMaterial()
    : UniaxialMaterial(0, MAT_TAG_RinaldinMasonryMaterial),
      modelType(Shear),
      parameterCount(9),
      elasticStiffness(0.0),
      yieldForce(0.0),
      firstPostYieldStiffness(0.0),
      maximumForce(0.0),
      secondPostYieldStiffness(0.0),
      ultimateDeformation(0.0),
      yieldDeformation(0.0),
      peakDeformation(0.0)
{
    std::fill(parameters, parameters + MaxParameterCount, 0.0);
    revertToStart();
}

const char *RinaldinMasonryMaterial::getClassType() const
{
    return modelType == Shear
        ? "RinaldinMasonryShear"
        : "RinaldinMasonryFlexural";
}

void RinaldinMasonryMaterial::setDerivedParameters()
{
    elasticStiffness = parameters[0];
    yieldForce = parameters[1];
    yieldDeformation = elasticStiffness > 0.0
        ? yieldForce / elasticStiffness
        : 0.0;

    firstPostYieldStiffness = parameters[2];
    maximumForce = parameters[3];
    secondPostYieldStiffness = modelType == Shear
        ? parameters[4]
        : parameters[9];
    ultimateDeformation = parameters[8];
    const double forcePeakDeformation = firstPostYieldStiffness > Tiny
        ? yieldDeformation
            + (maximumForce - yieldForce) / firstPostYieldStiffness
        : yieldDeformation;
    peakDeformation = modelType == Flexural
        ? std::min(forcePeakDeformation, ultimateDeformation)
        : forcePeakDeformation;
}

bool RinaldinMasonryMaterial::validateParameters(bool printMessage) const
{
    bool valid = true;
    auto report = [&](const char *message) {
        valid = false;
        if (printMessage) {
            opserr << "WARNING Rinaldin masonry material " << getTag()
                   << ": " << message << endln;
        }
    };

    if (parameterCount != (modelType == Shear ? 9 : 10)) {
        report("parameter count does not match the selected model");
    }
    for (int i = 0; i < parameterCount; ++i) {
        if (!std::isfinite(parameters[i])) {
            report("all parameters must be finite");
            break;
        }
    }
    if (elasticStiffness <= 0.0) {
        report("Kel must be greater than zero");
    }
    if (modelType == Shear) {
        if (yieldForce <= 0.0 || maximumForce < yieldForce) {
            report("require Fmax >= Fy > 0");
        }
        if (maximumForce > yieldForce && firstPostYieldStiffness <= 0.0) {
            report("K1pl must be positive when Fmax is greater than Fy");
        }
        if (ultimateDeformation <= 0.0 || peakDeformation > ultimateDeformation) {
            report("Uult must not be smaller than the deformation at Fmax");
        }
        const double CF = parameters[5];
        const double alpha = parameters[6];
        const double beta = parameters[7];
        if (CF < 0.0 || CF > 1.0) {
            report("CF must be in [0, 1]");
        }
        if (alpha <= 0.0 || alpha > 1.0) {
            report("alpha must be in (0, 1]");
        }
        if (beta < 0.0) {
            report("beta must be nonnegative");
        }
    } else {
        if (yieldForce <= 0.0 || maximumForce < yieldForce) {
            report("require Fmax >= Fy > 0");
        }
        if (maximumForce > yieldForce && firstPostYieldStiffness <= 0.0) {
            report("K1pl must be positive when Fmax is greater than Fy");
        }
        if (ultimateDeformation <= yieldDeformation) {
            report("Uult must be greater than the yield deformation Fy/Kel");
        }
        const double CC = parameters[4];
        const double CF = parameters[5];
        const double alpha = parameters[6];
        const double CD = parameters[7];
        if (CC <= 0.0) {
            report("CC must be greater than zero");
        }
        if (CF < 0.0 || CF > 1.0 || CD < 0.0) {
            report("CF must be in [0, 1] and CD must be nonnegative");
        }
        if (alpha <= 0.0 || alpha > 1.0) {
            report("alpha must be in (0, 1]");
        }
    }
    return valid;
}

double RinaldinMasonryMaterial::envelopeStress(double strain) const
{
    const double absoluteStrain = std::abs(strain);
    if (modelType == Shear && absoluteStrain > ultimateDeformation) {
        // Retain the residual shear strength beyond Uult instead of dropping to zero.
        const double residualForce = maximumForce
            + secondPostYieldStiffness
                * (ultimateDeformation - peakDeformation);
        return signOf(strain) * std::max(0.0, residualForce);
    }

    double force = 0.0;
    if (absoluteStrain <= yieldDeformation) {
        force = elasticStiffness * absoluteStrain;
    } else if (absoluteStrain <= peakDeformation) {
        force = yieldForce
            + firstPostYieldStiffness * (absoluteStrain - yieldDeformation);
    } else {
        const double transitionForce = modelType == Flexural
            ? yieldForce + firstPostYieldStiffness
                * (peakDeformation - yieldDeformation)
            : maximumForce;
        force = transitionForce
            + secondPostYieldStiffness * (absoluteStrain - peakDeformation);
        force = std::max(0.0, force);
    }
    return signOf(strain) * force;
}

double RinaldinMasonryMaterial::envelopeTangent(double strain) const
{
    const double absoluteStrain = std::abs(strain);
    if (modelType == Shear && absoluteStrain > ultimateDeformation) {
        return Tiny * elasticStiffness;
    }
    if (absoluteStrain < yieldDeformation) {
        return elasticStiffness;
    }
    if (absoluteStrain < peakDeformation) {
        return firstPostYieldStiffness;
    }
    return std::abs(secondPostYieldStiffness) > Tiny
        ? secondPostYieldStiffness
        : Tiny * elasticStiffness;
}

double RinaldinMasonryMaterial::unloadingStiffness(
    double maximumAbsoluteStrain) const
{
    if (modelType != Shear || maximumAbsoluteStrain <= yieldDeformation) {
        return elasticStiffness;
    }

    // Equations (18)-(21): degrade K1 linearly to alpha*K1 at Uult.
    const double alpha = parameters[6];
    const double denominator = ultimateDeformation - yieldDeformation;
    if (denominator <= Tiny) {
        return alpha * elasticStiffness;
    }
    const double ratio = std::clamp(
        (maximumAbsoluteStrain - yieldDeformation) / denominator,
        0.0,
        1.0);
    return elasticStiffness * (1.0 + (alpha - 1.0) * ratio);
}

double RinaldinMasonryMaterial::dissipatedEnergy(const State &state) const
{
    // Subtract secant strain energy from external work for Equation (22).
    return std::max(0.0, state.work - 0.5 * state.stress * state.strain);
}

bool RinaldinMasonryMaterial::reached(
    double value,
    double boundary,
    int direction) const
{
    return direction > 0 ? value >= boundary : value <= boundary;
}

double RinaldinMasonryMaterial::lineValue(
    double strain,
    double strain1,
    double stress1,
    double strain2,
    double stress2,
    double fallbackStiffness,
    double &tangent)
{
    const double delta = strain2 - strain1;
    tangent = std::abs(delta) > Tiny
        ? (stress2 - stress1) / delta
        : fallbackStiffness;
    return stress1 + tangent * (strain - strain1);
}

void RinaldinMasonryMaterial::startShearReversalPath(int newDirection)
{
    const double maximumAbsoluteStrain = std::max(
        std::abs(committed.maxPositiveStrain),
        std::abs(committed.minNegativeStrain));
    const double stiffness = unloadingStiffness(maximumAbsoluteStrain);
    const double CF = parameters[5];
    const double beta = parameters[7];

    trial.break1Stress = CF * committed.stress;
    trial.break1Strain = committed.strain
        + (trial.break1Stress - committed.stress) / stiffness;

    const double shift = maximumForce > Tiny
        ? beta * dissipatedEnergy(committed) / maximumForce
        : 0.0;
    trial.targetStrain = -committed.strain + newDirection * shift;
    // Connect the reversal path continuously to the opposite backbone.
    trial.targetStress = envelopeStress(trial.targetStrain);
    trial.break2Strain = trial.targetStrain;
    trial.break2Stress = trial.targetStress;
    trial.pathStage = 1;

    // Collapse to a stable line if a small cycle crosses the target.
    if (reached(trial.break1Strain, trial.targetStrain, newDirection)) {
        trial.break1Strain = trial.targetStrain;
        trial.break1Stress = trial.targetStress;
    }
}

void RinaldinMasonryMaterial::startFlexuralReversalPath(int newDirection)
{
    const double CC = parameters[4];
    const double CF = parameters[5];
    const double alpha = parameters[6];
    const double CD = parameters[7];
    const double oldSign = signOf(committed.stress);
    const double forceMagnitude = std::abs(committed.stress);

    // So.ph.i CF is the force drop relative to Fy, not the residual force ratio.
    const double firstStiffness = alpha * elasticStiffness;
    trial.break1Stress = oldSign * (
        forceMagnitude - CF * yieldForce);
    trial.break1Strain = committed.strain
        + (trial.break1Stress - committed.stress) / firstStiffness;

    // CD shifts the latest yield displacement to the right by CD*uy.
    const double secondStiffness = std::max(
        CC * firstPostYieldStiffness,
        Tiny * elasticStiffness);
    trial.break2Strain = oldSign * yieldDeformation
        + CD * yieldDeformation;
    trial.break2Stress = trial.break1Stress
        + secondStiffness * (
            trial.break2Strain - trial.break1Strain);

    // Branch 6 connects D to the elastic-limit point for reversed loading.
    trial.targetStrain = -oldSign * yieldDeformation;
    trial.targetStress = -oldSign * yieldForce;

    // Keep break points inside the path for roundoff or small cycles.
    if (reached(trial.break1Strain, trial.targetStrain, newDirection)) {
        trial.break1Strain = trial.targetStrain;
        trial.break1Stress = trial.targetStress;
        trial.break2Strain = trial.targetStrain;
        trial.break2Stress = trial.targetStress;
    } else if (reached(trial.break2Strain, trial.targetStrain, newDirection)) {
        trial.break2Strain = 0.5 * (
            trial.break1Strain + trial.targetStrain);
        trial.break2Stress = trial.break1Stress
            + secondStiffness * (trial.break2Strain - trial.break1Strain);
    }
    trial.pathStage = 1;
}

void RinaldinMasonryMaterial::startSecondaryReversalPath(int newDirection)
{
    // A reversal inside an existing unloading/reloading path returns to the
    // previously reached envelope point in the new loading direction. This
    // keeps stress continuous and avoids treating a Newton trial as a jump
    // directly back to the backbone.
    const double envelopeStrain = newDirection > 0
        ? committed.maxPositiveStrain
        : committed.minNegativeStrain;
    const bool envelopePointIsAhead = newDirection > 0
        ? envelopeStrain > committed.strain + Tiny
        : envelopeStrain < committed.strain - Tiny;

    if (!envelopePointIsAhead) {
        trial.pathStage = 0;
        return;
    }

    trial.break1Strain = envelopeStrain;
    trial.break1Stress = envelopeStress(envelopeStrain);
    trial.break2Strain = envelopeStrain;
    trial.break2Stress = trial.break1Stress;
    trial.targetStrain = envelopeStrain;
    trial.targetStress = trial.break1Stress;
    trial.pathStage = 1;
}

void RinaldinMasonryMaterial::startReversalPath(int newDirection)
{
    trial.direction = newDirection;
    trial.pathStartStrain = committed.strain;
    trial.pathStartStress = committed.stress;

    if (committed.pathStage != 0) {
        startSecondaryReversalPath(newDirection);
        return;
    }

    if (std::abs(committed.strain) <= yieldDeformation + Tiny) {
        trial.pathStage = 0;
        return;
    }

    if (modelType == Shear) {
        startShearReversalPath(newDirection);
    } else {
        startFlexuralReversalPath(newDirection);
    }
}

void RinaldinMasonryMaterial::evaluateCurrentPath(double strain)
{
    if (trial.pathStage == 0) {
        trial.stress = envelopeStress(strain);
        trial.tangent = envelopeTangent(strain);
        return;
    }

    if (!reached(strain, trial.break1Strain, trial.direction)) {
        trial.stress = lineValue(
            strain,
            trial.pathStartStrain,
            trial.pathStartStress,
            trial.break1Strain,
            trial.break1Stress,
            elasticStiffness,
            trial.tangent);
        return;
    }

    trial.pathStage = 2;
    if (!reached(strain, trial.break2Strain, trial.direction)) {
        trial.stress = lineValue(
            strain,
            trial.break1Strain,
            trial.break1Stress,
            trial.break2Strain,
            trial.break2Stress,
            Tiny * elasticStiffness,
            trial.tangent);
        return;
    }

    trial.pathStage = 3;
    if (!reached(strain, trial.targetStrain, trial.direction)) {
        trial.stress = lineValue(
            strain,
            trial.break2Strain,
            trial.break2Stress,
            trial.targetStrain,
            trial.targetStress,
            Tiny * elasticStiffness,
            trial.tangent);
        return;
    }

    trial.pathStage = 0;
    trial.stress = envelopeStress(strain);
    trial.tangent = envelopeTangent(strain);
}

int RinaldinMasonryMaterial::setTrialStrain(double strain, double strainRate)
{
    trial = committed;
    trial.strain = strain;
    const double increment = strain - committed.strain;
    if (std::abs(increment) <= Tiny) {
        return 0;
    }

    const int newDirection = increment > 0.0 ? 1 : -1;
    if (committed.direction == 0) {
        trial.direction = newDirection;
        trial.pathStage = 0;
    } else if (newDirection != committed.direction) {
        startReversalPath(newDirection);
    }

    evaluateCurrentPath(strain);
    trial.maxPositiveStrain = std::max(
        committed.maxPositiveStrain,
        strain);
    trial.minNegativeStrain = std::min(
        committed.minNegativeStrain,
        strain);
    trial.work = committed.work
        + 0.5 * (committed.stress + trial.stress) * increment;
    return 0;
}

double RinaldinMasonryMaterial::getStrain()
{
    return trial.strain;
}

double RinaldinMasonryMaterial::getStress()
{
    return trial.stress;
}

double RinaldinMasonryMaterial::getTangent()
{
    return trial.tangent;
}

double RinaldinMasonryMaterial::getInitialTangent()
{
    return elasticStiffness;
}

int RinaldinMasonryMaterial::commitState()
{
    committed = trial;
    return 0;
}

int RinaldinMasonryMaterial::revertToLastCommit()
{
    trial = committed;
    return 0;
}

int RinaldinMasonryMaterial::revertToStart()
{
    State initial{};
    initial.tangent = elasticStiffness > 0.0 ? elasticStiffness : 0.0;
    committed = initial;
    trial = initial;
    return 0;
}

UniaxialMaterial *RinaldinMasonryMaterial::getCopy()
{
    auto *copy = new RinaldinMasonryMaterial(
        getTag(), modelType, parameters, parameterCount);
    copy->committed = committed;
    copy->trial = trial;
    return copy;
}

int RinaldinMasonryMaterial::sendSelf(
    int commitTag,
    Channel &theChannel)
{
    Vector data(29);
    data(0) = getTag();
    data(1) = static_cast<int>(modelType);
    data(2) = parameterCount;
    for (int i = 0; i < MaxParameterCount; ++i) {
        data(3 + i) = parameters[i];
    }
    data(13) = committed.strain;
    data(14) = committed.stress;
    data(15) = committed.tangent;
    data(16) = committed.maxPositiveStrain;
    data(17) = committed.minNegativeStrain;
    data(18) = committed.work;
    data(19) = committed.direction;
    data(20) = committed.pathStage;
    data(21) = committed.pathStartStrain;
    data(22) = committed.pathStartStress;
    data(23) = committed.break1Strain;
    data(24) = committed.break1Stress;
    data(25) = committed.break2Strain;
    data(26) = committed.break2Stress;
    data(27) = committed.targetStrain;
    data(28) = committed.targetStress;

    if (getDbTag() == 0) {
        setDbTag(theChannel.getDbTag());
    }
    return theChannel.sendVector(getDbTag(), commitTag, data);
}

int RinaldinMasonryMaterial::recvSelf(
    int commitTag,
    Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    Vector data(29);
    const int result = theChannel.recvVector(getDbTag(), commitTag, data);
    if (result < 0) {
        opserr << "RinaldinMasonryMaterial::recvSelf failed" << endln;
        return result;
    }

    setTag(static_cast<int>(data(0)));
    modelType = static_cast<ModelType>(static_cast<int>(data(1)));
    parameterCount = static_cast<int>(data(2));
    for (int i = 0; i < MaxParameterCount; ++i) {
        parameters[i] = data(3 + i);
    }
    setDerivedParameters();
    committed.strain = data(13);
    committed.stress = data(14);
    committed.tangent = data(15);
    committed.maxPositiveStrain = data(16);
    committed.minNegativeStrain = data(17);
    committed.work = data(18);
    committed.direction = static_cast<int>(data(19));
    committed.pathStage = static_cast<int>(data(20));
    committed.pathStartStrain = data(21);
    committed.pathStartStress = data(22);
    committed.break1Strain = data(23);
    committed.break1Stress = data(24);
    committed.break2Strain = data(25);
    committed.break2Stress = data(26);
    committed.targetStrain = data(27);
    committed.targetStress = data(28);
    trial = committed;
    return 0;
}

void RinaldinMasonryMaterial::Print(OPS_Stream &stream, int flag)
{
    stream << getClassType() << ", tag: " << getTag() << endln;
    stream << "  Kel: " << elasticStiffness
           << ", Fy: " << yieldForce
           << ", K1pl: " << firstPostYieldStiffness
           << ", Fmax: " << maximumForce
           << ", K2pl: " << secondPostYieldStiffness
           << ", Uult: " << ultimateDeformation << endln;
}
