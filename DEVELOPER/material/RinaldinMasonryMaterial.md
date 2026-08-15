# Rinaldin masonry uniaxial materials

本模块将 Rinaldin 等（2016）第 2.4 节的砌体剪切与弯曲滞回规则实现为两个
OpenSees `UniaxialMaterial`。材料只描述一个相对自由度上的力-变形关系，需与
`zeroLength` 单元配合使用。

## 剪切材料

```tcl
uniaxialMaterial RinaldinMasonryShear $tag \
    $Kel $Fy $K1pl $Fmax $K2pl $CF $alpha $beta $Uult
```

- `Kel`：弹性刚度。
- `Fy`：开裂/屈服力。
- `K1pl`：第一后弹性分支刚度。
- `Fmax`：第一后弹性分支末端的最大力。
- `K2pl`：第二后弹性分支刚度；软化时为负值。
- `CF`：从骨架卸载时的力比参数。
- `alpha`：极限变形处卸载刚度与弹性刚度之比。
- `beta`：基于滞回耗能的强度退化系数。
- `Uult`：极限位移。

## 弯曲材料

```tcl
uniaxialMaterial RinaldinMasonryFlexural $tag \
    $Kel $Fy $K1pl $Fmax $CC $CF $alpha $CD $Uult $K2pl
```

- `Kel`：弹性转动刚度。
- `Fy`：开裂弯矩。
- `K1pl`：第一后弹性转动刚度。
- `Fmax`：最大弯矩。
- `CC`：第二卸载分支相对后弹性分支的刚度系数。
- `CF`：离开骨架后第一卸载分支的力降比例。
- `alpha`：第一卸载分支相对弹性分支的刚度系数。
- `CD`：D 点相对于最新屈服位移的向右移动比例，`u_D=sign(M_A)*(1+CD)*Fy/Kel`。
- `Uult`：极限转角。
- `K2pl`：达到 `Fmax` 或 `Uult` 中较早者后采用的第二后弹性分支刚度；论文双折线模型通常取 `0`。

## zeroLength 示例

```tcl
# 中部剪切弹簧沿局部平移自由度 1 工作。
uniaxialMaterial RinaldinMasonryShear 1 \
    1.6 18.6 0.3 25.0 -1.2 0.8 0.8 0.06 34.0
element zeroLength 101 11 12 -mat 1 -dir 1

# 端部弯曲弹簧绕转动自由度 5 工作。
uniaxialMaterial RinaldinMasonryFlexural 2 \
    1.6 18.6 0.3 25.0 1.0 0.1 0.8 0.1 34.0 0.0
element zeroLength 102 21 22 -mat 2 -dir 5
```

两个端部弯曲弹簧应使用两个材料标签，使其分别保存各自的循环历史。轴向弹簧可
单独使用 `Elastic` 材料。模型必须通过约束明确处理未工作的自由度。

## 当前验证范围

- MinGW GCC 14 与 MSVC 2022 均已通过源文件编译。
- 已提供 Tcl 回归脚本，用于检查骨架关键点和循环过程中的有限数值；需在完整的新建 `OpenSees.exe` 中运行。
- `CF`、`CC`、`alpha`、`CD` 的实现依据论文图 5-7 与 So.ph.i. 2021 手册。
- 在获得 So.ph.i. 同位移历程的逐点输出后，还应追加数值对照回归测试。
