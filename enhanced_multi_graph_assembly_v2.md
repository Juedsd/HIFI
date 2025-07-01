# 多图整合基因组组装系统设计文档 (架构重构版)

## 1. 项目概述

### 1.1 项目背景
基于现有 hifiasm 基因组组装器，构建一个支持多数据源和多图结构整合的增强型组装系统。通过整合 HiFi reads、Ultra-long ONT reads 和参考基因组信息，结合智能化的GAP填充策略，提高组装质量和连续性。

### 1.2 项目目标
- 构建 HiFi reads 双向 String Graph（主图）
- 支持参考基因组辅助图（HiFi→RefGraph）
- 支持 Ultra-long ONT 辅助图（ONT→HiFi ULGraph, ONT→RefGraph）
- 实现智能化的多图合并策略
- 智能GAP填充模块，支持多种填充策略
- **架构重构：模块化设计，符合hifiasm现有代码结构**

## 2. 架构重构设计

### 2.1 重构后的模块分层架构

```
┌─────────────────────────────────────────────────────────────┐
│                    Assembly.cpp (主工作流集成)                │
├─────────────────────────────────────────────────────────────┤
│ • perform_multi_graph_mapping()                            │
│ • multi_graph_assembly_pipeline()                          │
│ • 集成各模块功能，控制整体流程                                │
└─────────────────┬───────────────────────────────────────────┘
                  │
    ┌─────────────┼─────────────┐
    │             │             │
┌───▼────┐  ┌────▼─────┐  ┌────▼─────┐
│inter.cpp│  │gchain_map │  │Overlaps  │
│(比对算法)│  │.cpp(图合并)│  │.cpp(数据) │
├─────────┤  ├──────────┤  ├──────────┤
│HiFi→Ref │  │图权重管理 │  │4图overlap │
│ONT→Ref  │  │图结构优化 │  │统计管理   │
│复用UL算法│  │复用gchain │  │质量验证   │
└─────────┘  └──────────┘  └──────────┘
                  │
            ┌─────▼─────┐
            │ref_genome │
            │.cpp(基础) │
            ├───────────┤
            │数据管理   │
            │FASTA加载  │
            │UL索引构建 │
            └───────────┘
```

### 2.2 模块职责重新划分

#### **ref_genome.cpp** - 基础数据管理层
**核心职责**：参考基因组数据结构管理
```cpp
// ✅ 保留在此模块的功能
├── ref_genome_init/destroy                    // 数据结构管理
├── ref_genome_load_fasta                      // FASTA文件加载
├── ref_genome_convert_to_all_reads           // 转换为All_reads格式
├── ref_genome_build_ul_index                 // UL索引构建
├── ref_genome_save/load_cache                // 缓存管理
└── ref_config_t 相关配置函数                  // 配置管理
```

#### **inter.cpp** - 比对算法扩展层
**核心职责**：复用UL比对算法处理多类型比对
```cpp
// 🆕 新增到此模块的功能
├── hifi_map_to_reference()                   // HiFi→Ref比对 (复用ul_map_lchain)
├── ont_map_to_reference()                    // ONT→Ref比对 (复用UL算法)
├── minimizer_search_in_reference()           // minimizer查找 (复用现有)
├── validate_reference_alignment()            // 比对质量验证
└── 复用现有UL比对核心算法                     // mg_lchain_gen, ul_map_lchain等
```

#### **gchain_map.cpp** - 图算法扩展层
**核心职责**：复用图链接算法处理多图合并
```cpp
// 🆕 新增到此模块的功能
├── merge_main_and_ref_graphs()               // 图合并 (复用mg_gchain1_dp)
├── apply_multi_graph_weights()               // 多图权重管理
├── optimize_ref_graph_structure()            // 图结构优化
├── resolve_graph_conflicts()                 // 图冲突解决
└── 复用现有图链接核心算法                     // mg_gchain1_dp, 动态规划等
```

#### **Overlaps.cpp** - 数据管理扩展层
**核心职责**：管理多类型overlap数据
```cpp
// 🆕 新增到此模块的功能
├── four_graph_overlaps_t 管理函数            // init/destroy/统计
├── print_mapping_statistics()                // 多图统计输出
├── validate_multi_graph_quality()            // 多图质量验证
├── export_multi_graph_results()              // 结果导出
└── gap_fill_* 相关函数                      // GAP填充算法
```

#### **Assembly.cpp** - 主工作流集成层
**核心职责**：集成各模块，控制整体工作流
```cpp
// 🆕 新增到此模块的功能
├── perform_multi_graph_mapping()             // 主工作流控制
├── multi_graph_assembly_pipeline()           // 完整组装管道
├── 调用各模块的算法函数                       // 模块间协调
└── 错误处理和资源管理                        // 异常处理
```

## 3. 核心数据结构设计

### 3.1 参考基因组结构 (ref_genome.cpp)
```cpp
typedef struct {
    char *fasta_path;                    // FASTA文件路径
    uint32_t n_seq;                      // 序列数量
    ref_chromosome_t *chromosomes;       // 染色体数组
    char *merged_seq;                    // 合并序列
    uint64_t total_length;               // 总长度
    uint64_t n_bases;                    // 碱基数
    uint64_t n_chunks;                   // 分块数
    
    All_reads *ref_reads;                // 转换为All_reads格式
    ul_idx_t *ul_index;                  // UL索引
    int index_built;                     // 索引构建状态
} ref_genome_t;
```

### 3.2 四图overlap管理结构 (Overlaps.cpp)
```cpp
typedef struct {
    ma_hit_t_alloc *hifi_hifi_overlaps;  // HiFi↔HiFi主图
    ma_hit_t_alloc *ont_hifi_overlaps;   // ONT→HiFi辅助图
    ma_hit_t_alloc *hifi_ref_overlaps;   // HiFi→Ref辅助图
    ma_hit_t_alloc *ont_ref_overlaps;    // ONT→Ref辅助图
    
    uint64_t n_hifi_hifi_edges;          // 各类边的数量统计
    uint64_t n_ont_hifi_edges;
    uint64_t n_hifi_ref_edges;
    uint64_t n_ont_ref_edges;
    
    multi_graph_weights_t weights;        // 权重配置
} four_graph_overlaps_t;
```

### 3.3 多图权重配置 (gchain_map.cpp)
```cpp
typedef struct {
    double hifi_hifi_weight;             // 主图权重 (默认1.0)
    double ont_hifi_weight;              // ONT辅助权重 (默认0.3)
    double hifi_ref_weight;              // 参考辅助权重 (默认0.5)
    double ont_ref_weight;               // ONT参考权重 (默认0.2)
    double identity_bonus;               // 相似度奖励
    double length_penalty;               // 长度惩罚
} multi_graph_weights_t;
```

## 4. 算法复用策略

### 4.1 UL比对算法复用 (inter.cpp)

#### HiFi→Ref比对实现
```cpp
int hifi_map_to_reference(const All_reads *hifi_reads,
                         ref_genome_t *ref,
                         ma_hit_t_alloc *ref_overlaps,
                         const hifiasm_opt_t *opt) {
    
    // 复用现有的UL比对核心算法
    for (uint64_t i = 0; i < hifi_reads->total_reads; i++) {
        
        // 1. 复用 ul_map_lchain 进行链接
        ul_map_lchain(ab, i, Get_READ(*hifi_reads, i), 
                     Get_READ_LENGTH(*hifi_reads, i), 
                     opt->mz_win, opt->k_mer_length,
                     ref->ul_index, overlap_list, cl, 
                     opt->bw_thres, opt->max_n_chain, ...);
        
        // 2. 转换UL overlap到标准overlap格式
        convert_ul_overlaps_to_ma_hits(ul_overlaps, ref_overlaps);
    }
    
    return 0;
}
```

#### 复用的核心函数映射
| 新功能 | 复用的hifiasm函数 | 位置 |
|--------|------------------|------|
| HiFi→Ref minimizer查找 | `mg_lchain_gen` | inter.cpp |
| 链接扩展算法 | `ul_map_lchain` | gchain_map.cpp |
| 动态规划链接 | `mg_gchain1_dp` | inter.cpp |
| 质量评估 | `mg_gchain_set_mapq` | inter.cpp |

### 4.2 图合并算法复用 (gchain_map.cpp)

#### 多图合并实现
```cpp
int merge_main_and_ref_graphs(asg_t *main_graph,
                             ma_hit_t_alloc *ref_overlaps,
                             const multi_graph_weights_t *weights,
                             asg_t **merged_graph) {
    
    // 1. 复制主图作为基础
    *merged_graph = asg_dup(main_graph);
    
    // 2. 基于参考overlaps添加约束边
    for (uint64_t i = 0; i < ref_overlaps->length; i++) {
        ma_hit_t *hit = &ref_overlaps->buffer[i];
        
        // 3. 复用gchain质量评估
        if (evaluate_overlap_with_gchain_quality(hit, *merged_graph)) {
            // 4. 添加权重调整后的边
            add_weighted_edge(*merged_graph, hit, weights);
        }
    }
    
    // 5. 复用现有图清理算法
    asg_cleanup(*merged_graph);
    
    return 0;
}
```

## 5. GAP填充模块设计

### 5.1 GAP填充集成位置
**集成到 Overlaps.cpp**，在unitig construction之后执行：

```cpp
// Overlaps.cpp 中添加GAP填充功能
typedef enum {
    GAP_FILL_NONE = 0,           // 不填充
    GAP_FILL_REF_SEQ,            // 使用参考序列填充
    GAP_FILL_LOWERCASE,          // 使用小写序列填充
    GAP_FILL_N_PADDING           // 使用N字符填充
} gap_fill_mode_t;

typedef struct {
    uint64_t gap_start;          // GAP起始位置
    uint64_t gap_end;            // GAP结束位置
    uint64_t gap_size;           // GAP大小
    gap_fill_mode_t fill_mode;   // 填充模式
    float confidence;            // 填充置信度
    char *fill_sequence;         // 填充序列
    int validation_status;       // 验证状态
} gap_info_t;
```

### 5.2 GAP填充算法
```cpp
int perform_intelligent_gap_filling(asg_t *graph,
                                   ref_genome_t *ref_genome,
                                   const hifiasm_opt_t *opt) {
    
    // 1. 检测GAP
    gap_info_t *gaps = detect_gaps_in_graph(graph);
    
    // 2. 分类GAP并选择填充策略
    for (each gap) {
        if (gap_size <= opt->min_gap_size) continue;
        if (gap_size >= opt->max_gap_size) continue;
        
        // 智能选择填充模式
        gap->fill_mode = select_gap_fill_strategy(gap, ref_genome, opt);
    }
    
    // 3. 执行填充
    execute_gap_filling(gaps, graph, ref_genome);
    
    // 4. 质量验证
    validate_gap_filling_quality(gaps, graph);
    
    return 0;
}
```

## 6. 命令行接口

### 6.1 新增参数 (CommandLines.cpp)
```bash
# 基本多图选项
--ref FILE              参考基因组FASTA文件
--ul FILE               Ultra-long ONT reads文件
--enable-multi-graph    启用多图模式 (默认关闭)

# 权重控制
--hifi-ref-weight FLOAT HiFi→Ref图权重 (默认0.5)
--ont-hifi-weight FLOAT ONT→HiFi图权重 (默认0.3)
--ont-ref-weight FLOAT  ONT→Ref图权重 (默认0.2)

# GAP填充
--gap-fill MODE         GAP填充模式: none|ref|lowercase|n (默认none)
--min-gap INT           最小GAP大小 (默认10)
--max-gap INT           最大GAP大小 (默认10000)
--gap-threshold FLOAT   GAP填充置信度阈值 (默认0.7)

# 输出控制
--verbose-multi         详细多图信息输出
--multi-prefix STR      多图输出文件前缀
```

### 6.2 使用示例
```bash
# 1. 基本多图组装
./hifiasm -o output --ref ref.fasta --enable-multi-graph hifi.fq

# 2. 三图组装 (HiFi + ONT + Ref)
./hifiasm -o output --ref ref.fasta --ul ont.fq --enable-multi-graph hifi.fq

# 3. 带GAP填充的组装
./hifiasm -o output --ref ref.fasta --gap-fill ref --gap-threshold 0.8 hifi.fq

# 4. 自定义权重组装
./hifiasm -o output --ref ref.fasta --ul ont.fq \
    --hifi-ref-weight 0.6 --ont-hifi-weight 0.4 \
    --verbose-multi hifi.fq
```

## 7. 开发阶段规划

### 7.1 第一阶段：模块重构 (2周)
**目标**：重新组织现有代码，建立正确的模块架构

**主要任务**：
- 从 `ref_genome.cpp` 迁移比对算法到 `inter.cpp`
- 从 `ref_genome.cpp` 迁移图合并算法到 `gchain_map.cpp`
- 从 `ref_genome.cpp` 迁移overlap管理到 `Overlaps.cpp`
- 在 `Assembly.cpp` 中建立主工作流

**文件修改**：
```cpp
// inter.cpp - 新增比对函数
+ int hifi_map_to_reference(...)
+ int ont_map_to_reference(...)

// gchain_map.cpp - 新增图合并函数  
+ int merge_main_and_ref_graphs(...)
+ int apply_multi_graph_weights(...)

// Overlaps.cpp - 新增overlap管理
+ four_graph_overlaps_t* init_four_graph_overlaps(...)
+ void print_mapping_statistics(...)

// Assembly.cpp - 新增主工作流
+ int perform_multi_graph_mapping(...)

// ref_genome.cpp - 精简为数据管理
- 移除比对和图合并相关函数
+ 保留数据结构管理函数
```

### 7.2 第二阶段：算法实现 (3周)
**目标**：实现核心比对和图合并算法

**主要任务**：
- 在 `inter.cpp` 中实现基于UL算法的HiFi→Ref比对
- 在 `gchain_map.cpp` 中实现基于现有算法的图合并
- 在 `Overlaps.cpp` 中完善overlap数据管理
- 实现质量验证和统计功能

### 7.3 第三阶段：GAP填充模块 (2周)
**目标**：在 `Overlaps.cpp` 中实现GAP填充功能

**主要任务**：
- GAP检测算法
- 智能填充策略选择
- 多种填充模式实现
- 填充质量验证

### 7.4 第四阶段：集成测试 (2周)
**目标**：完整工作流测试和优化

**主要任务**：
- 模块间集成测试
- 性能优化
- 错误处理完善
- 文档和使用指南

## 8. 预期效果

### 8.1 架构优势
- **模块化清晰**：每个模块职责单一，易于维护
- **代码复用最大化**：充分利用hifiasm现有算法
- **扩展性好**：便于后续功能扩展
- **符合hifiasm风格**：与现有代码结构一致

### 8.2 性能目标
- **组装连续性提升**：scaffold N50 提升 ≥ 20%
- **计算效率**：总运行时间增长 ≤ 25%
- **内存使用**：内存占用增长 ≤ 30%
- **GAP填充准确率**：≥ 95%

### 8.3 功能完整性
- ✅ 支持三种图构建模式 (HiFi, HiFi+Ref, HiFi+ONT+Ref)
- ✅ 智能权重调整和冲突解决
- ✅ 多种GAP填充策略
- ✅ 详细的质量评估和统计
- ✅ 灵活的配置选项

## 9. 风险评估与应对

### 9.1 技术风险
**风险**：模块重构可能引入新的bug
**应对**：分阶段迁移，每个阶段充分测试

**风险**：算法复用可能不完全兼容
**应对**：深入理解现有算法，必要时进行适配

### 9.2 性能风险
**风险**：多图处理增加计算复杂度
**应对**：充分利用现有并行框架，优化关键路径

## 10. 总结

通过重新设计模块架构，将功能正确分配到各个模块中，既保持了与hifiasm现有代码结构的一致性，又实现了多图整合的功能目标。这种设计方案具有以下优势：

1. **架构合理**：符合hifiasm的模块化设计理念
2. **复用充分**：最大化利用现有算法和数据结构
3. **扩展性强**：为未来功能扩展奠定良好基础
4. **维护性好**：模块职责清晰，便于长期维护

该设计既实现了多图整合基因组组装的技术目标，又确保了代码的工程质量和可维护性。