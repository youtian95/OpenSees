// ============================================================================
// 文件: MasonryBendingMat.cpp
// 功能: MasonryBendingMat 砌体弯曲单轴材料实现 + OPS 工厂函数
// 说明: 样板代码(工厂函数/状态管理/并行通讯)已完成,
//       核心本构规律在 setTrialStrain 中以 TODO 标出, 由开发者手动实现
// 用法: 编译后在脚本中使用:
//       Tcl:        uniaxialMaterial MasonryBendingMat matTag? Ke? Mmax? <可选参数>
//       OpenSeesPy: ops.uniaxialMaterial('MasonryBendingMat', matTag, Ke, Mmax, ...)
// 输入: matTag(整数); Ke Mmax(必选实数); uu~gamma2(可选, 尾部省略用默认值)
// 输出: 创建 MasonryBendingMat 材料对象并注册到模型构建器
// ============================================================================

#include "MasonryBendingMat.h"

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>

// 已创建的材料实例数(用于只在首次创建时打印署名横幅)
static int numMasonryBendingMat = 0;

// 工厂函数: uniaxialMaterial 命令的入口
// 命令语法(< > 内为可选参数, 尾部省略时使用构造函数的默认参数):
//   uniaxialMaterial MasonryBendingMat matTag? Ke? Mmax? <uu?> <R_My?> <alpha1?> <alpha2?> <gamma1?> <gamma2?>
void *OPS_MasonryBendingMat(void)
{
  // 首次创建时打印一次署名信息
  if (numMasonryBendingMat == 0) {
    opserr << "MasonryBendingMat uniaxial material - Written by youtian95\n";
    numMasonryBendingMat = 1;
  }

  // 至少需要 matTag + Ke + Mmax 三个参数
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING invalid args, want: uniaxialMaterial MasonryBendingMat tag? Ke? Mmax? <uu?> <R_My?> <alpha1?> <alpha2?> <gamma1?> <gamma2?>\n";
    return 0;
  }

  // 1) 读取整数参数: matTag
  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MasonryBendingMat tag" << endln;
    return 0;
  }

  // 2) 读取实数参数: 只读用户实际给出的个数(封顶8个),
  //    再按实际个数调用构造函数, 尾部未给出的参数由头文件中的默认参数补齐
  double dData[8];
  int numArgs = OPS_GetNumRemainingInputArgs();
  numData = (numArgs > 8) ? 8 : numArgs;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid MasonryBendingMat parameters\n";
    return 0;
  }

  // 3) 按实际读到的参数个数创建材料, 缺省参数自动生效
  UniaxialMaterial *theMaterial = 0;
  switch (numData) {
    case 2:  theMaterial = new MasonryBendingMat(iData[0], dData[0], dData[1]); break;
    case 3:  theMaterial = new MasonryBendingMat(iData[0], dData[0], dData[1], dData[2]); break;
    case 4:  theMaterial = new MasonryBendingMat(iData[0], dData[0], dData[1], dData[2], dData[3]); break;
    case 5:  theMaterial = new MasonryBendingMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]); break;
    case 6:  theMaterial = new MasonryBendingMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]); break;
    case 7:  theMaterial = new MasonryBendingMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]); break;
    default: theMaterial = new MasonryBendingMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7]); break;
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MasonryBendingMat\n";
    return 0;
  }
  return theMaterial;
}

// 全参构造函数
// tag: 材料标签; 其余参数含义见头文件
// 基类必须用 (tag, MAT_TAG_MasonryBendingMat) 初始化
MasonryBendingMat::MasonryBendingMat(int tag, double _Ke, double _Mmax,
                                     double _uu, double _R_My,
                                     double _alpha1, double _alpha2,
                                     double _gamma1, double _gamma2)
  : UniaxialMaterial(tag, MAT_TAG_MasonryBendingMat),
    Ke(_Ke), Mmax(_Mmax), uu(_uu), R_My(_R_My),
    alpha1(_alpha1), alpha2(_alpha2), gamma1(_gamma1), gamma2(_gamma2),
    trialStrain(0.0), trialStrainRate(0.0),
    trialStress(0.0), trialTangent(_Ke),
    cStrain(0.0), cStress(0.0), cTangent(_Ke)
{
  // TODO: 如有其他派生参数/历史变量, 在此初始化
}

// 默认构造函数: 供 FEM_ObjectBroker 在并行/数据库恢复时使用
MasonryBendingMat::MasonryBendingMat()
  : UniaxialMaterial(0, MAT_TAG_MasonryBendingMat),
    Ke(0.0), Mmax(0.0), uu(0.008), R_My(0.7),
    alpha1(0.2), alpha2(0.1), gamma1(1.2), gamma2(1.2),
    trialStrain(0.0), trialStrainRate(0.0),
    trialStress(0.0), trialTangent(0.0),
    cStrain(0.0), cStress(0.0), cTangent(0.0)
{
}

// 析构函数
MasonryBendingMat::~MasonryBendingMat()
{
}

// 设置试应变并计算应力与切线刚度
// strain: 试应变; strainRate: 试应变率
// ★★★ 核心本构规律写在这里: 由 trialStrain 与历史变量计算 trialStress/trialTangent ★★★
int
MasonryBendingMat::setTrialStrain(double strain, double strainRate)
{
  trialStrain = strain;
  trialStrainRate = strainRate;

  // TODO: 实现弯曲本构规律, 例如:
  //   - 判断加载/卸载方向
  //   - 计算骨架曲线应力(弯矩)
  //   - 计算切线刚度
  //   - 更新历史变量
  trialStress  = 0.0;  // TODO: 替换为本构计算的应力(弯矩)
  trialTangent = Ke;   // TODO: 替换为本构计算的切线刚度

  return 0;
}

// 以下四个 get 函数返回当前试状态
double MasonryBendingMat::getStrain(void)         { return trialStrain; }
double MasonryBendingMat::getStress(void)         { return trialStress; }
double MasonryBendingMat::getTangent(void)        { return trialTangent; }
double MasonryBendingMat::getInitialTangent(void) { return Ke; }

// 提交状态: 试状态转为已提交状态
int
MasonryBendingMat::commitState(void)
{
  cStrain  = trialStrain;
  cStress  = trialStress;
  cTangent = trialTangent;
  // TODO: 一并提交你的历史变量
  return 0;
}

// 回滚到上次提交状态(迭代不收敛时调用)
int
MasonryBendingMat::revertToLastCommit(void)
{
  trialStrain  = cStrain;
  trialStress  = cStress;
  trialTangent = cTangent;
  // TODO: 一并回滚你的历史变量
  return 0;
}

// 回滚到初始状态(重新开始分析时调用)
int
MasonryBendingMat::revertToStart(void)
{
  trialStrain = 0.0;
  trialStrainRate = 0.0;
  trialStress = 0.0;
  trialTangent = Ke;
  cStrain = 0.0;
  cStress = 0.0;
  cTangent = Ke;
  // TODO: 一并重置你的历史变量
  return 0;
}

// 返回材料完整副本(单元迭代时调用)
// 成员均为值类型, 默认拷贝构造即可; 若后续加入指针成员需自行深拷贝
UniaxialMaterial *
MasonryBendingMat::getCopy(void)
{
  return new MasonryBendingMat(*this);
}

// 打包发送(并行/数据库): 材料tag + 8个材料参数 + 已提交状态
// commitTag: 提交标签; theChannel: 通讯通道
int
MasonryBendingMat::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(12);
  data(0) = this->getTag();
  data(1) = Ke;      data(2) = Mmax;
  data(3) = uu;      data(4) = R_My;
  data(5) = alpha1;  data(6) = alpha2;
  data(7) = gamma1;  data(8) = gamma2;
  data(9) = cStrain;  data(10) = cStress;  data(11) = cTangent;
  // TODO: 若增加了历史变量, 需扩大 data 并一并打包

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0)
    opserr << "MasonryBendingMat::sendSelf() - failed to send data\n";
  return res;
}

// 接收恢复(并行/数据库): 打包顺序必须与 sendSelf 严格对应
int
MasonryBendingMat::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(12);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "MasonryBendingMat::recvSelf() - failed to receive data\n";
    return res;
  }
  this->setTag((int)data(0));
  Ke = data(1);      Mmax = data(2);
  uu = data(3);      R_My = data(4);
  alpha1 = data(5);  alpha2 = data(6);
  gamma1 = data(7);  gamma2 = data(8);
  cStrain = data(9);  cStress = data(10);  cTangent = data(11);
  // 用已提交状态初始化试状态
  trialStrain = cStrain;
  trialStress = cStress;
  trialTangent = cTangent;
  return 0;
}

// 打印材料信息
// s: 输出流; flag: 输出详细程度标志
void
MasonryBendingMat::Print(OPS_Stream &s, int flag)
{
  s << "MasonryBendingMat tag: " << this->getTag() << endln;
  s << "  Ke: " << Ke << " Mmax: " << Mmax << " uu: " << uu << endln;
  s << "  R_My: " << R_My << " alpha1: " << alpha1 << " alpha2: " << alpha2
    << " gamma1: " << gamma1 << " gamma2: " << gamma2 << endln;
  // TODO: 按需打印更多状态信息
}
