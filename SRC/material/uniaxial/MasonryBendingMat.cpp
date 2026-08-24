// ============================================================================
// 文件: MasonryBendingMat.cpp
// 功能: MasonryBendingMat 砌体弯曲单轴材料实现与 OPS 工厂函数
// 使用流程: 工厂函数读取材料参数 -> 状态机判断骨架/卸载/再加载分支
//           -> 计算弯矩和切线 -> commitState 提交历史
// 输入: matTag, Ke, Mmax, uu 以及可选的 R_My、CF、CD、gamma1、gamma2
// 输出: OpenSees UniaxialMaterial，返回弯矩-转角关系
// ============================================================================

#include "MasonryBendingMat.h"

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>

#include <algorithm>
#include <cmath>

// 已创建的材料实例数，用于只在首次创建时打印署名信息。
static int numMasonryBendingMat = 0;

// 工厂函数：读取命令参数并创建 MasonryBendingMat。
void *OPS_MasonryBendingMat(void)
{
  if (numMasonryBendingMat == 0) {
    opserr << "MasonryBendingMat uniaxial material - Written by youtian95\n";
    numMasonryBendingMat = 1;
  }

  // tag + Ke + Mmax 至少需要三个输入参数，uu 及其后的参数可省略。
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING invalid args, want: uniaxialMaterial MasonryBendingMat "
           << "tag? Ke? Mmax? <uu?> <R_My?> <CF?> <CD?> "
           << "<gamma1?> <gamma2?>\n";
    return 0;
  }

  int tag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid MasonryBendingMat tag" << endln;
    return 0;
  }

  // 两个必选实数加六个可选实数，共八个实数参数。
  double data[8];
  numData = OPS_GetNumRemainingInputArgs();
  if (numData > 8)
    numData = 8;
  if (OPS_GetDoubleInput(&numData, data) != 0) {
    opserr << "WARNING invalid MasonryBendingMat parameters" << endln;
    return 0;
  }

  UniaxialMaterial *material = 0;
  switch (numData) {
    case 2: material = new MasonryBendingMat(tag, data[0], data[1]); break;
    case 3: material = new MasonryBendingMat(tag, data[0], data[1], data[2]); break;
    case 4: material = new MasonryBendingMat(tag, data[0], data[1], data[2], data[3]); break;
    case 5: material = new MasonryBendingMat(tag, data[0], data[1], data[2], data[3], data[4]); break;
    case 6: material = new MasonryBendingMat(tag, data[0], data[1], data[2], data[3], data[4], data[5]); break;
    case 7: material = new MasonryBendingMat(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6]); break;
    default: material = new MasonryBendingMat(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]); break;
  }

  if (material == 0) {
    opserr << "WARNING could not create MasonryBendingMat" << endln;
    return 0;
  }
  return material;
}

// 全参构造函数。
MasonryBendingMat::MasonryBendingMat(
    int tag, double _Ke, double _Mmax, double _uu, double _R_My,
    double _CF, double _CD, double _gamma1, double _gamma2)
  : UniaxialMaterial(tag, MAT_TAG_MasonryBendingMat),
    Ke(_Ke), initialMmax(_Mmax), Mmax(_Mmax), committedMmax(_Mmax), uu(_uu), R_My(_R_My),
    CF(_CF), CD(_CD), gamma1(_gamma1), gamma2(_gamma2),
    tState(), cState(), My(0.0), uy(0.0), Kp(0.0)
{
  updateDerivedParameters();
  tState = initialState();
  cState = initialState();
}

// 默认构造函数：供 FEM_ObjectBroker 并行/数据库恢复使用。
MasonryBendingMat::MasonryBendingMat()
  : UniaxialMaterial(0, MAT_TAG_MasonryBendingMat),
    Ke(0.0), initialMmax(0.0), Mmax(0.0), committedMmax(0.0), uu(0.008), R_My(0.7),
    CF(0.2), CD(0.1), gamma1(1.2), gamma2(1.2),
    tState(), cState(), My(0.0), uy(0.0), Kp(0.0)
{
  tState = initialState();
  cState = initialState();
}

MasonryBendingMat::~MasonryBendingMat()
{
}

// 设置试转角，并按“先判断分支、后执行规则”的顺序计算试状态。
int MasonryBendingMat::setTrialStrain(double strain, double strainRate)
{
  tState = cState;
  tState.strain = strain;
  tState.strainRate = strainRate;

  // 根据当前试转角相对已提交转角的增量确定加载方向。
  const double du = tState.strain - cState.strain;
  if (du > 1.0e-14)
    tState.ldir = 1.0;
  else if (du < -1.0e-14)
    tState.ldir = -1.0;
  else
    tState.ldir = cState.ldir;

  // 状态转换必须先于规则函数执行，规则函数依赖 tState.branch。
  checkTransitions();

  // cState.branch 表示从哪个已提交分支出发，tState.branch 表示下一分支。
  switch (cState.branch) {
    case ELASTIC: ruleFromElastic(); break;
    case HARDENING_1: ruleFromHardening1(); break;
    case UNLOADING_1: ruleFromUnloading1(); break;
    case UNLOADING_2: ruleFromUnloading2(); break;
    case UNLOADING_3: ruleFromUnloading3(); break;
    case RELOADING_FROM_UNLOADING: ruleFromReloading(); break;
  }

  return 0;
}

double MasonryBendingMat::getStrain(void) { return tState.strain; }
double MasonryBendingMat::getStress(void) { return tState.stress; }
double MasonryBendingMat::getTangent(void) { return tState.tangent; }
double MasonryBendingMat::getInitialTangent(void) { return Ke; }

// 设置当前试算轴力对应的弯曲骨架峰值。
// Mmax: 当前试算最大弯矩；应在本次 setTrialStrain() 之前调用。
int MasonryBendingMat::setTrialBackbone(double trialMmax)
{
  // 用初始承载力的百万分之一作为数值下限，避免零承载力导致派生参数除零。
  const double minimumMmax = 1.0e-6 * initialMmax;
  Mmax = std::max(trialMmax, minimumMmax);
  updateDerivedParameters();
  return 0;
}

// 提交试状态。
int MasonryBendingMat::commitState(void)
{
  cState = tState;
  committedMmax = Mmax;
  return 0;
}

// 回滚到上一收敛状态。
int MasonryBendingMat::revertToLastCommit(void)
{
  tState = cState;
  Mmax = committedMmax;
  updateDerivedParameters();
  return 0;
}

// 回滚到初始状态。
int MasonryBendingMat::revertToStart(void)
{
  Mmax = initialMmax;
  committedMmax = initialMmax;
  updateDerivedParameters();
  tState = initialState();
  cState = initialState();
  return 0;
}

// 返回完整材料副本。
UniaxialMaterial *MasonryBendingMat::getCopy(void)
{
  return new MasonryBendingMat(*this);
}

// 发送材料参数和完整已提交状态。
int MasonryBendingMat::sendSelf(int commitTag, Channel &theChannel)
{
  // data(0)为材料tag，data(1..8)为参数，data(9)为已提交骨架峰值，data(10..18)为状态。
  static Vector data(19);
  data(0) = getTag();
  data(1) = Ke; data(2) = initialMmax; data(3) = uu; data(4) = R_My;
  data(5) = CF; data(6) = CD; data(7) = gamma1; data(8) = gamma2;
  data(9) = committedMmax;
  data(10) = cState.strain;
  data(11) = cState.strainRate;
  data(12) = cState.stress;
  data(13) = cState.tangent;
  data(14) = static_cast<int>(cState.branch);
  data(15) = cState.ldir;
  data(16) = cState.revStrain;
  data(17) = cState.revStress;
  data(18) = cState.directUnloading3 ? 1.0 : 0.0;

  int result = theChannel.sendVector(getDbTag(), commitTag, data);
  if (result < 0)
    opserr << "MasonryBendingMat::sendSelf() - failed to send data\n";
  return result;
}

// 接收材料参数和完整已提交状态。
int MasonryBendingMat::recvSelf(
    int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(19);
  int result = theChannel.recvVector(getDbTag(), commitTag, data);
  if (result < 0) {
    opserr << "MasonryBendingMat::recvSelf() - failed to receive data\n";
    return result;
  }

  setTag(static_cast<int>(data(0)));
  Ke = data(1); initialMmax = data(2); uu = data(3); R_My = data(4);
  CF = data(5); CD = data(6); gamma1 = data(7); gamma2 = data(8);
  committedMmax = data(9);
  Mmax = committedMmax;
  updateDerivedParameters();

  cState.strain = data(10);
  cState.strainRate = data(11);
  cState.stress = data(12);
  cState.tangent = data(13);
  cState.branch = static_cast<Branch>(static_cast<int>(data(14)));
  cState.ldir = data(15);
  cState.revStrain = data(16);
  cState.revStress = data(17);
  cState.directUnloading3 = data(18) != 0.0;
  tState = cState;
  return 0;
}

// 打印材料信息。
void MasonryBendingMat::Print(OPS_Stream &s, int flag)
{
  s << "MasonryBendingMat tag: " << getTag() << endln;
  s << "  Ke: " << Ke << " initialMmax: " << initialMmax << " Mmax: " << Mmax << " uu: " << uu << endln;
  s << "  My: " << My << " uy: " << uy << " Kp: " << Kp << endln;
  s << "  R_My: " << R_My << " CF: " << CF
    << " CD: " << CD << " gamma1: " << gamma1
    << " gamma2: " << gamma2 << endln;
}

// 根据当前试算最大弯矩更新所有依赖骨架峰值的派生参数。
void MasonryBendingMat::updateDerivedParameters()
{
  My = R_My * Mmax;
  uy = Ke != 0.0 ? My / Ke : 0.0;
  Kp = uu > uy ? (Mmax - My) / (uu - uy) : 0.0;
}

// 计算无退化的双折线弯矩-转角骨架。
void MasonryBendingMat::backbone(double rotation, double &moment, double &tangent)
{
  const double absoluteRotation = std::abs(rotation);
  const double sign = rotation < 0.0 ? -1.0 : 1.0;

  if (absoluteRotation <= uy) {
    moment = Ke * rotation;
    tangent = Ke;
  } else if (absoluteRotation <= uu) {
    moment = sign * (My + Kp * (absoluteRotation - uy));
    tangent = Kp;
  } else {
    // 达到极限转角后保持最后的弯矩水平，不引入额外退化。
    moment = sign * Mmax;
    tangent = 0.0;
  }
}

// 根据反转点A临时计算三段卸载的两个折点，不写入历史状态。
// state: 当前试状态; unload1: CF确定的第一折点; unload2: CD确定的B点
void MasonryBendingMat::computeUnloadingPoints(const State &state, double &unload1Strain, double &unload1Stress, double &unload2Strain, double &unload2Stress)
{
  if (state.directUnloading3) {
    unload1Strain = state.revStrain;
    unload1Stress = state.revStress;
    unload2Strain = state.revStrain;
    unload2Stress = state.revStress;
    return;
  }

  const double k1 = gamma1 * Ke;
  const double k2 = gamma2 * Kp;
  const double safeK1 = std::abs(k1) > 1.0e-14 ? k1 : Ke;
  const double safeK2 = std::abs(k2) > 1.0e-14 ? k2 : Kp;

  // CF控制第一折点的弯矩下降比例。
  unload1Stress = (1.0 - CF) * state.revStress;
  unload1Strain = state.revStrain + (unload1Stress - state.revStress) / safeK1;

  // CD控制第二折点B相对屈服转角的水平位置，B点弯矩由第二段直线连续计算。
  const double side = state.revStress >= 0.0 ? 1.0 : -1.0;
  unload2Strain = side * (1.0 + CD) * uy;
  unload2Stress = unload1Stress + safeK2 * (unload2Strain - unload1Strain);
}

// 根据反转点A的弯矩方向计算第三段卸载目标点。
// 正弯矩反转时冲向负向屈服点，负弯矩反转时冲向正向屈服点。
double MasonryBendingMat::unloadingTargetStrain(const State &state) const
{
  // 从再加载直接进入第三段时，目标点沿新的加载方向确定。
  if (state.directUnloading3)
    return state.ldir > 0.0 ? uy : -uy;

  // 普通卸载进入第三段时，目标点取反转点弯矩的反向屈服点。
  if (state.revStress > 1.0e-14)
    return -uy;
  if (state.revStress < -1.0e-14)
    return uy;

  // 反转点弯矩接近零时，使用当前加载方向作为退化情形下的后备判断。
  return state.ldir > 0.0 ? uy : -uy;
}

// 按A-B-C-反向屈服点的三段折线计算卸载路径。
// state: 当前试状态，rotation: 当前试转角
void MasonryBendingMat::evaluateUnloadingPath(State &state, double rotation)
{
  const double k1 = gamma1 * Ke;
  const double k2 = gamma2 * Kp;
  const double safeK1 = std::abs(k1) > 1.0e-14 ? k1 : Ke;
  const double safeK2 = std::abs(k2) > 1.0e-14 ? k2 : Kp;
  double unload1Strain = 0.0;
  double unload1Stress = 0.0;
  double unload2Strain = 0.0;
  double unload2Stress = 0.0;
  computeUnloadingPoints(state, unload1Strain, unload1Stress, unload2Strain, unload2Stress);

  if (state.branch == UNLOADING_1) {
    state.stress = state.revStress + safeK1 * (rotation - state.revStrain);
    state.tangent = safeK1;
  } else if (state.branch == UNLOADING_2) {
    state.stress = unload1Stress + safeK2 * (rotation - unload1Strain);
    state.tangent = safeK2;
  } else if (state.branch == UNLOADING_3) {
    const double targetStrain = unloadingTargetStrain(state);
    double targetStress = 0.0;
    double targetTangent = 0.0;
    backbone(targetStrain, targetStress, targetTangent);
    const double denominator = targetStrain - unload2Strain;
    state.tangent = std::abs(denominator) > 1.0e-14 ? (targetStress - unload2Stress) / denominator : safeK2;
    state.stress = unload2Stress + state.tangent * (rotation - unload2Strain);
  }
}

// 判断当前试状态是否跨越了下一分支边界。
bool MasonryBendingMat::checkTransitions()
{
  const double du = tState.strain - cState.strain;
  if (std::abs(du) <= 1.0e-14)
    return false;

  if (cState.branch == ELASTIC) {
    if (std::abs(tState.strain) > uy) {
      tState.branch = HARDENING_1;
      return true;
    }
    return false;
  }

  if (cState.branch == HARDENING_1) {
    if (tState.ldir * cState.strain < 0.0) {
      tState.branch = UNLOADING_1;
      return true;
    }
    return false;
  }

  if (cState.branch == UNLOADING_1) {
    // 卸载段内只要增量方向反转，就进入反向再加载。
    // 不能用当前弯矩正负判断，否则再加载尚未过零时再次反向会被误判。
    if (tState.ldir * cState.ldir < 0.0) {
      tState.branch = RELOADING_FROM_UNLOADING;
      return true;
    }

    const double assumedStress = cState.stress + cState.tangent * du;
    const double direction = tState.ldir > 0.0 ? 1.0 : -1.0;
    const double unload1Stress = (1.0 - CF) * cState.revStress;
    const double side = cState.revStress >= 0.0 ? 1.0 : -1.0;
    const double unload2Strain = side * (1.0 + CD) * uy;

    // 只要当前位移已经越过B点，就直接进入第三段，不再经过第二段。
    if ((tState.strain - unload2Strain) * direction > 0.0) {
      tState.branch = UNLOADING_3;
      return true;
    }

    if ((assumedStress - unload1Stress) * direction > 0.0) {
      tState.branch = UNLOADING_2;
      return true;
    }
    return false;
  }

  if (cState.branch == UNLOADING_2) {
    if (tState.ldir * cState.ldir < 0.0) {
      tState.branch = RELOADING_FROM_UNLOADING;
      return true;
    }
    const double direction = tState.ldir > 0.0 ? 1.0 : -1.0;
    const double side = cState.revStress >= 0.0 ? 1.0 : -1.0;
    const double unload2Strain = side * (1.0 + CD) * uy;
    if ((tState.strain - unload2Strain) * direction > 0.0) {
      tState.branch = UNLOADING_3;
      return true;
    }
    return false;
  }

  if (cState.branch == UNLOADING_3) {
    if (tState.ldir * cState.ldir < 0.0) {
      tState.branch = RELOADING_FROM_UNLOADING;
      return true;
    }
    const double targetStrain = unloadingTargetStrain(cState);
    if ((tState.strain - targetStrain) * tState.ldir > 0.0) {
      tState.branch = HARDENING_1;
      return true;
    }
    return false;
  }

  if (cState.branch == RELOADING_FROM_UNLOADING) {
    // 必须先判断再次反向，否则反向增量可能被误判为越过再加载目标点。
    if (tState.ldir * cState.ldir < 0.0) {
      // 只有反转点已经超过上一加载侧的B点，才重新经历三段卸载；
      // 尚未达到B点时，直接连接到反向屈服点。
      const double bStrain = cState.ldir * (1.0 + CD) * uy;
      if ((cState.strain - bStrain) * cState.ldir > 0.0) {
        tState.branch = UNLOADING_1;
      } else {
        tState.branch = UNLOADING_3;
      }
      return true;
    }
    // 再加载沿K1直线上升，在当前加载方向与骨架相交时回到骨架。
    const double K1 = gamma1 * Ke;
    const double assumedStress = cState.stress + K1 * du;
    double backboneStress = 0.0;
    double backboneTangent = 0.0;
    backbone(tState.strain, backboneStress, backboneTangent);
    if ((assumedStress - backboneStress) * cState.ldir >= 0.0) {
      tState.branch = HARDENING_1;
      return true;
    }
    return false;
  }

  return false;
}

// 从弹性分支进入弹性或屈服后骨架。
void MasonryBendingMat::ruleFromElastic()
{
  if (tState.branch == ELASTIC || tState.branch == HARDENING_1)
    backbone(tState.strain, tState.stress, tState.tangent);
}

// 从屈服后骨架进入骨架或第一段卸载。
void MasonryBendingMat::ruleFromHardening1()
{
  if (tState.branch == HARDENING_1) {
    backbone(tState.strain, tState.stress, tState.tangent);
  } else if (tState.branch == UNLOADING_1) {
    tState.revStrain = cState.strain;
    tState.revStress = cState.stress;
    tState.directUnloading3 = false;
    tState.tangent = gamma1 * Ke;
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  }
}

// 从第一段卸载进入第一段、第二段或反向再加载。
void MasonryBendingMat::ruleFromUnloading1()
{
  if (tState.branch == UNLOADING_1) {
    tState.tangent = cState.tangent;
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  } else if (tState.branch == UNLOADING_2 ||
             tState.branch == UNLOADING_3) {
    evaluateUnloadingPath(tState, tState.strain);
  } else if (tState.branch == RELOADING_FROM_UNLOADING) {
    tState.tangent = gamma1 * Ke;
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  }
}

// 从第二段卸载进入第二段、第三段、反向再加载或骨架。
void MasonryBendingMat::ruleFromUnloading2()
{
  if (tState.branch == UNLOADING_2 || tState.branch == UNLOADING_3) {
    evaluateUnloadingPath(tState, tState.strain);
  } else if (tState.branch == RELOADING_FROM_UNLOADING) {
    tState.tangent = gamma1 * Ke;
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  }
}

// 从第三段卸载进入第三段、反向再加载或反向骨架。
void MasonryBendingMat::ruleFromUnloading3()
{
  if (tState.branch == UNLOADING_3) {
    evaluateUnloadingPath(tState, tState.strain);
  } else if (tState.branch == RELOADING_FROM_UNLOADING) {
    tState.tangent = gamma1 * Ke;
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  } else if (tState.branch == HARDENING_1) {
    backbone(tState.strain, tState.stress, tState.tangent);
  }
}

// 从反向再加载进入再加载、再次卸载或骨架。
void MasonryBendingMat::ruleFromReloading()
{
  if (tState.branch == RELOADING_FROM_UNLOADING) {
    tState.tangent = gamma1 * Ke;
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  } else if (tState.branch == HARDENING_1) {
    backbone(tState.strain, tState.stress, tState.tangent);
  } else if (tState.branch == UNLOADING_1) {
    tState.revStrain = cState.strain;
    tState.revStress = cState.stress;
    tState.directUnloading3 = false;
    tState.tangent = gamma1 * Ke;
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  } else if (tState.branch == UNLOADING_3) {
    // 再加载在零力之前反向：以当前点作为第三段起点，直接连接反向屈服点。
    tState.revStrain = cState.strain;
    tState.revStress = cState.stress;
    tState.directUnloading3 = true;
    const double targetStrain = unloadingTargetStrain(tState);
    double targetStress = 0.0;
    double targetTangent = 0.0;
    backbone(targetStrain, targetStress, targetTangent);
    const double denominator = targetStrain - tState.revStrain;
    tState.tangent = (targetStress - tState.revStress) / denominator;
    tState.stress = tState.revStress + tState.tangent * (tState.strain - tState.revStrain);
  }
}

// 返回初始状态。
MasonryBendingMat::State MasonryBendingMat::initialState() const
{
  State state = {};
  state.tangent = Ke;
  return state;
}
