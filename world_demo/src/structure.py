"""
核心数据结构定义
WorldNode 和 GenomePointer
"""

import uuid
from typing import List, Any, Optional, Union


# 常量配置
R_SW = 5
R_MW = 15
R_LW = 30
SEQ_LEN = 250
SKELETON_SAMPLE_SIZE = 1000000


class BioSequence:
    """基础数据单元，封装原始字符串"""
    
    def __init__(self, seq_id: str, seq: str):
        """
        Args:
            seq_id: 唯一标识符
            seq: DNA 序列字符串 (A, T, C, G)
        """
        self.id = seq_id
        self.seq = seq
    
    def __repr__(self):
        """返回序列摘要"""
        seq_preview = self.seq[:20] + "..." if len(self.seq) > 20 else self.seq
        return f"BioSequence(id={self.id}, seq={seq_preview}, len={len(self.seq)})"


class GenomePointer:
    """基因组指针（内存模拟版本）
    用于存储序列坐标，模拟从 2-bit 格式提取序列
    """
    
    def __init__(self, chrom: str, pos: int, seq: str):
        """
        Args:
            chrom: 染色体标识
            pos: 位置坐标
            seq: 序列内容（内存中，模拟从 2-bit 提取）
        """
        self.chrom = chrom
        self.pos = pos
        self.seq = seq  # 内存中的序列引用
    
    def get_sequence(self) -> str:
        """获取序列内容"""
        return self.seq
    
    def __repr__(self):
        return f"GenomePointer(chrom={self.chrom}, pos={self.pos}, len={len(self.seq)})"


class WorldNode:
    """DAG 索引节点（v3.0）
    支持 Multi-Parent 结构，使用 UUID 进行去重
    """
    
    def __init__(self, center_ptr: Union[BioSequence, GenomePointer], radius: int, layer_level: int):
        """
        Args:
            center_ptr: 中心序列指针（BioSequence 或 GenomePointer）
            radius: 固定半径
            layer_level: 层级标识 (1=SW, 2=MW, 3=LW)
        """
        # 生成唯一 UUID 用于搜索去重
        layer_name = {1: "SW", 2: "MW", 3: "LW"}.get(layer_level, "UNK")
        self.node_id = f"{layer_name}_{uuid.uuid4().hex[:8]}"
        
        self.center_ptr = center_ptr
        self.radius = radius
        self.layer = layer_level
        
        # 子节点引用（支持 Multi-Parent）
        self.children: List[WorldNode] = []
        
        # 叶子节点数据（仅 SW 有效）
        self.disk_offset = -1  # 磁盘偏移（当前为内存模拟，暂不使用）
        self.data_count = 0  # 叶子节点数据数量
    
    def get_center_sequence(self) -> str:
        """获取中心序列字符串"""
        if isinstance(self.center_ptr, BioSequence):
            return self.center_ptr.seq
        elif isinstance(self.center_ptr, GenomePointer):
            return self.center_ptr.get_sequence()
        else:
            return str(self.center_ptr)
    
    def add_child(self, child: 'WorldNode'):
        """添加子节点（支持 Multi-Parent）"""
        self.children.append(child)
    
    def is_overlapping(self, query_seq: Union[BioSequence, str], query_radius: int, distance_func) -> bool:
        """
        判断查询球体与节点球体是否相交
        
        Args:
            query_seq: Query 序列对象或字符串
            query_radius: Query 容错半径
            distance_func: 距离计算函数
            
        Returns:
            如果 dist(query, self.center) <= self.radius + query_radius，返回 True
        """
        # 获取查询序列字符串
        if isinstance(query_seq, BioSequence):
            query_str = query_seq.seq
        elif isinstance(query_seq, GenomePointer):
            query_str = query_seq.get_sequence()
        else:
            query_str = str(query_seq)
        
        # 获取中心序列字符串
        center_str = self.get_center_sequence()
        
        # 创建临时 BioSequence 对象用于距离计算
        query_obj = BioSequence("temp_query", query_str)
        center_obj = BioSequence("temp_center", center_str)
        
        dist = distance_func(query_obj, center_obj)
        return dist <= self.radius + query_radius
    
    def __repr__(self):
        layer_name = {1: "SW", 2: "MW", 3: "LW"}.get(self.layer, "UNK")
        return f"WorldNode(id={self.node_id}, layer={layer_name}, radius={self.radius}, children={len(self.children)}, data_count={self.data_count})"
