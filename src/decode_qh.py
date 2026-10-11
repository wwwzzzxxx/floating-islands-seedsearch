#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
decode_qh.py — 从 Cubiomes Viewer 的预计算表解码出全部「四联女巫小屋」基

背景：四联小屋的基是一个静态可预计算的集合，Cubiomes Viewer 直接把它当数据文件
随软件发布（仓库里的 rc/qh/，28 个文件、差分编码、约 4.1 MB）。
本仓库初版是用 CUDA 穷举整个 2^48 空间算出来的（RTX 4090 上 740 秒），
后来发现可以直接读表——已交叉验证：解码结果 755,745 个，与 GPU 穷举数量完全一致，
且本仓库 1483 个种子的小屋基 1483/1483 全部命中。

用法：
    python decode_qh.py <cubiomes-viewer 的 rc/qh 目录> [-o 输出文件] [--salt 14357620]

    # 拿到表（只需 rc/qh 目录，不用克隆整个仓库）：
    git clone --depth 1 --filter=blob:none --sparse https://github.com/Cubitect/cubiomes-viewer.git
    cd cubiomes-viewer && git sparse-checkout set rc/qh && cd ..

    python decode_qh.py cubiomes-viewer/rc/qh -o hutbases.txt

输出格式：每行一个 48 位基（十进制），已按数值排序。
--salt 默认 14357620（MC 1.18.2 的 Swamp_Hut salt）；换 MC 版本时用
    getStructureConfig(Swamp_Hut, <你的版本>, &sc) 取到的 sc.salt。
"""

import argparse
import io
import os
import sys

DEFAULT_SALT = 14357620  # MC_1_18_2 的 Swamp_Hut salt


def decode_qh(qh_dir, salt):
    """解码 rc/qh/ 目录，返回 48 位基的集合（已转换为原始种子空间）。"""
    bases = []
    for fn in sorted(os.listdir(qh_dir)):
        try:
            low = int(fn, 16)          # 文件名 = 低 20 位（salt 空间）
        except ValueError:
            continue                    # 跳过非十六进制文件名的杂物
        path = os.path.join(qh_dir, fn)
        if not os.path.isfile(path):
            continue
        mid = 0
        with io.open(path, "r", encoding="ascii", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    diff = int(line, 16)
                except ValueError:
                    continue
                if diff == 0:
                    break               # 0 是结束标记
                mid += diff             # 差分累加得高 28 位
                s48 = (mid << 20) + low
                bases.append((s48 - salt) & ((1 << 48) - 1))
    return bases


def main():
    ap = argparse.ArgumentParser(
        description="从 Cubiomes Viewer 的 rc/qh/ 表解码全部四联女巫小屋基")
    ap.add_argument("qh_dir", help="cubiomes-viewer 仓库的 rc/qh 目录")
    ap.add_argument("-o", "--out", default="hutbases.txt", help="输出文件（默认 hutbases.txt）")
    ap.add_argument("--salt", type=int, default=DEFAULT_SALT,
                    help="Swamp_Hut 的 salt，默认 %(default)s（MC 1.18.2）")
    args = ap.parse_args()

    if not os.path.isdir(args.qh_dir):
        sys.exit("目录不存在：%s" % args.qh_dir)

    bases = decode_qh(args.qh_dir, args.salt)
    if not bases:
        sys.exit("没解出任何基，请确认 %s 是 cubiomes-viewer 的 rc/qh 目录" % args.qh_dir)

    uniq = sorted(set(bases))
    with io.open(args.out, "w", encoding="ascii") as f:
        for b in uniq:
            f.write("%d\n" % b)

    print("解码完成：%d 个基（去重后 %d 个）-> %s" % (len(bases), len(uniq), args.out))
    print("提示：这些基的「小屋四联」锚在原点区域；要用到别处就用 "
          "moveStructure(base, regX, regZ)。")


if __name__ == "__main__":
    main()
