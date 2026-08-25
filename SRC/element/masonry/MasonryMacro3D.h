// 模块功能：定义三维局部平面内砌体宏单元 MasonryMacro3D 的 OpenSees 接口与数据结构。
// 使用流程：通过 element masonryMacro 命令创建单元，框架负责节点、材料、坐标变换、状态和并行通讯。
// 输入输出：输入两个六自由度节点、三类单轴材料和坐标变换；输出 12x12 刚度与 12 维恢复力。

#ifndef MasonryMacro3D_h
#define MasonryMacro3D_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class CrdTransf;
class Domain;
class FEM_ObjectBroker;
class Information;
class Node;
class OPS_Stream;
class Renderer;
class Response;
class UniaxialMaterial;

// 三维两节点砌体宏单元：每个节点具有六个自由度，但单元只在局部 x-y 平面内工作。
class MasonryMacro3D : public Element
{
public:
    // 创建供分析使用的三维局部平面内砌体宏单元。
    // tag: 单元标签。
    // nodeI: I 端节点标签。
    // nodeJ: J 端节点标签。
    // bendingMaterial: 弯曲材料原型，单元为 I、J 两端分别复制一份。
    // shearMaterial: 剪切材料原型。
    // axialMaterial: 轴向材料原型。
    // coordinateTransformation: 三维坐标变换原型，局部 x-y 为单元工作的平面。
    // width: pier 墙面内宽度。
    // thickness: pier 墙厚。
    // compressiveStrength: 砌体抗压强度。
    // cohesion: 砌体灰缝黏聚力。
    // diagonalTensileStrength: 砌体对角抗拉强度。
    // contraflexureDistance: pier 端部到反弯点的距离。
    // stressBlockCoefficient: 矩形压应力块系数。
    // frictionCoefficient: 灰缝摩擦系数。
    // maximumIterations: 内部剪切变形局部 Newton 迭代的最大次数。
    // relativeTolerance: 内部平衡残差的相对收敛容差。
    MasonryMacro3D(int tag, int nodeI, int nodeJ, UniaxialMaterial &bendingMaterial, UniaxialMaterial &shearMaterial, UniaxialMaterial &axialMaterial, CrdTransf &coordinateTransformation, double width = 0.0, double thickness = 0.0, double compressiveStrength = 0.0, double cohesion = 0.0, double diagonalTensileStrength = 0.0, double contraflexureDistance = 0.0, double stressBlockCoefficient = 0.85, double frictionCoefficient = 0.4, int maximumIterations = 30, double relativeTolerance = 1.0e-10);

    // 创建供 ObjectBroker 接收数据使用的空对象。
    MasonryMacro3D();

    // 释放单元拥有的材料副本和坐标变换副本。
    ~MasonryMacro3D();

    const char *getClassType(void) const { return "MasonryMacro3D"; }

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);

    // theDomain: 单元所属的 Domain；传入空指针表示解除关联。
    void setDomain(Domain *theDomain);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Vector &getResistingForce(void);

    // commitTag: 当前提交标签。
    // theChannel: 接收序列化数据的通道。
    int sendSelf(int commitTag, Channel &theChannel);

    // commitTag: 当前提交标签。
    // theChannel: 提供序列化数据的通道。
    // theBroker: 根据 classTag 创建材料和坐标变换的对象工厂。
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // theViewer: 图形渲染器。
    // displayMode: 显示模式。
    // fact: 位移放大系数。
    // modes: 可选显示模式名称数组。
    // numModes: modes 数组长度。
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes = 0, int numModes = 0);

    // s: 输出流。
    // flag: 输出格式标志。
    void Print(OPS_Stream &s, int flag = 0);

    // argv: 响应请求参数。
    // argc: 响应请求参数数量。
    // output: 响应描述输出流。
    Response *setResponse(const char **argv, int argc, OPS_Stream &output);

    // responseID: setResponse 分配的响应编号。
    // information: 承载响应结果的对象。
    int getResponse(int responseID, Information &information);

private:
    // 将全局矩阵投影到局部 x-y 平面的活动自由度，删除局部 z 向平动、绕 x 扭转和绕 y 转动分量。
    // matrix: 需要原位投影的 12x12 全局矩阵。
    int projectToActivePlane(Matrix &matrix) const;

    // 将全局向量投影到局部 x-y 平面的活动自由度。
    // vector: 需要原位投影的 12 维全局向量。
    int projectToActivePlane(Vector &vector) const;

    // 根据当前轴向压力计算弯曲和剪切峰值，并更新三个横向材料的试算骨架。
    // axialForce: 轴向材料当前试算力，拉伸为正、压缩为负。
    int updateMaterialBackbones(double axialForce);

    // 先用局部 Newton 法求解满足弯矩平衡的剪切弹簧变形，失败后改用局部括区间二分法，并把收敛变形写入三个横向材料的 trial 状态。
    // thetaI: I 端相对于单元弦的基本试算转角，逆时针为正。
    // thetaJ: J 端相对于单元弦的基本试算转角，逆时针为正。
    // axialForce: 当前轴力，拉伸为正。
    int solveInternalShearDeformation(double thetaI, double thetaJ, double axialForce);

    // 计算指定剪切变形下的局部平衡残差，并更新三个横向材料的 trial 状态。
    // thetaI: I 端相对于单元弦的基本试算转角，逆时针为正。
    // thetaJ: J 端相对于单元弦的基本试算转角，逆时针为正。
    // axialForce: 当前轴力，拉伸为正。
    // shearDeformation: 剪切弹簧试算变形。
    // residual: 返回局部平衡残差。
    // residualScale: 返回用于相对收敛判断的残差量级。
    int evaluateInternalShearResidual(double thetaI, double thetaJ, double axialForce, double shearDeformation, double &residual, double &residualScale);

    // 使用普通 Newton 法求解内部剪切变形。
    // thetaI: I 端相对于单元弦的基本试算转角，逆时针为正。
    // thetaJ: J 端相对于单元弦的基本试算转角，逆时针为正。
    // axialForce: 当前轴力，拉伸为正。
    // initialShearDeformation: Newton 迭代的初始剪切变形。
    int solveInternalShearDeformationByNewton(double thetaI, double thetaJ, double axialForce, double initialShearDeformation);

    // 从初始剪切变形附近寻找异号区间，并使用二分法求解内部剪切变形。
    // thetaI: I 端相对于单元弦的基本试算转角，逆时针为正。
    // thetaJ: J 端相对于单元弦的基本试算转角，逆时针为正。
    // axialForce: 当前轴力，拉伸为正。
    // initialShearDeformation: 二分搜索的中心剪切变形。
    int solveInternalShearDeformationByBisection(double thetaI, double thetaJ, double axialForce, double initialShearDeformation);

    ID connectedExternalNodes;
    Node *theNodes[2];
    CrdTransf *theCoordTransf;

    UniaxialMaterial *bendingMaterials[2];
    UniaxialMaterial *shearMaterial;
    UniaxialMaterial *axialMaterial;

    // 轴力相关承载力计算参数；width 为零时不启用轴力相关骨架更新。
    double width;
    double thickness;
    double compressiveStrength;
    double cohesion;
    double diagonalTensileStrength;
    double contraflexureDistance;
    double stressBlockCoefficient;
    double frictionCoefficient;

    // 内部剪切变形局部求解参数。
    int maximumIterations;
    double relativeTolerance;

    // 以下对象作为 OpenSees const 引用返回值的持久缓存。
    Matrix tangentStiffness;
    Matrix initialStiffness;
    Vector resistingForce;
};

#endif
