#!/usr/bin/env python3
"""
arg_list: generate all combinations (Cartesian product) from lists using a template.

Usage:
    arg_list [options] "<template>"

Placeholders:
   %i%                     - 1-based output line number
   %Nf%                    - whole item from list N (N≥1)
   %NfX M%                 - field M from list N using separator X (e.g., %1f,2%)
   %NfX1M1X2M2...%         - chained splits (e.g., %2f|1;2%)

Options:
   -f N        first N output lines
   -l N        last N output lines
   -L a,b      lines a to b (inclusive, 1-based)

Example:
    arg_list "%1f%: %2f|1%"
"""

import sys
import re
import argparse
from itertools import product
from typing import List, Iterable

# ---------------------------- Core logic ----------------------------

PLACEHOLDER_RE = re.compile(r"%((?P<file>\d+)f(?P<ops>(?:[^\d]\d+)*)?)%|%i%")

from typing import List


def generate_number_strings(start: int = 1, end: int = 10, step: int = 1) -> List[str]:
    """
    生成从 start 到 end（包含）的等差整数序列的字符串列表，步长为 step。

    参数:
        start (int): 起始数字（默认 1）
        end (int): 结束数字（默认 10）
        step (int): 步长，必须为正整数（默认 1）

    返回:
        List[str]: 字符串形式的数字列表

    示例:
        generate_number_strings(1, 5, 1)   → ["1", "2", "3", "4", "5"]
        generate_number_strings(0, 10, 2)  → ["0", "2", "4", "6", "8", "10"]
        generate_number_strings(5, 1, 1)   → []  （start > end 且 step > 0）
        generate_number_strings(10, 1, -1) → 报错（不支持负步长）
    """
    if step <= 0:
        raise ValueError("step must be a positive integer")
    if start > end:
        return []
    return [str(i) for i in range(start, end + 1, step)]


def write_lines_to_file(lines: List[str], filepath: str, mode: str = "a") -> None:
    """
    将字符串列表写入文本文件，每行一个元素。

    参数:
        lines (List[str]): 要写入的字符串列表。
        filepath (str): 目标文件路径。
        mode (str): 写入模式，默认 'a'（追加，不覆写）；
                    设为 'w' 可强制覆盖文件。

    示例:
        write_lines_to_file(["new entry"], "log.txt")          # 追加（默认）
        write_lines_to_file(["fresh start"], "data.txt", 'w')  # 覆盖
    """
    if not isinstance(lines, list):
        raise TypeError("lines must be a list of strings")
    if not all(isinstance(line, str) for line in lines):
        raise ValueError("All items in lines must be strings")

    try:
        with open(filepath, mode, encoding="utf-8") as f:
            for line in lines:
                f.write(line + "\n")
    except OSError as e:
        raise RuntimeError(f"Failed to write to file '{filepath}': {e}") from e


def read_lines_from_file(filepath: str) -> List[str]:
    """
    从文本文件中读取所有行，返回字符串列表（每行末尾的换行符已被去除）。

    参数:
        filepath (str): 要读取的文件路径。

    返回:
        List[str]: 文件中的每一行作为一个字符串元素。

    异常:
        FileNotFoundError: 文件不存在。
        OSError: 其他文件读取错误（如权限问题）。

    示例:
        lines = read_lines_from_file("data.txt")
        # 若 data.txt 内容为：
        #   apple
        #   banana
        # 则 lines == ["apple", "banana"]
    """
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n\r") for line in f]
        return lines
    except OSError as e:
        raise RuntimeError(f"Failed to read file '{filepath}': {e}") from e


def apply_ops(content: str, ops: str) -> str:
    """Apply chained split operations like ;1|2 on content."""
    if not ops:
        return content
    steps = re.findall(r"([^\d])(\d+)", ops)
    for sep, idx_str in steps:
        parts = content.split(sep)
        idx = int(idx_str) - 1
        if 0 <= idx < len(parts):
            content = parts[idx]
        else:
            return ""
    return content


def substitute(template: str, combo: List[str], line_no: int) -> str:
    """Expand template for one combination."""

    def repl(m: re.Match):
        if m.group(0) == "%i%":
            return str(line_no)
        file_no = int(m.group("file")) - 1
        if not (0 <= file_no < len(combo)):
            return ""
        content = combo[file_no]
        return apply_ops(content, m.group("ops") or "")

    return PLACEHOLDER_RE.sub(repl, template)


def select_indices(
    total: int, first: int = None, last: int = None, line_range: str = None
) -> set:
    """Return set of 1-based indices to output based on CLI options."""
    if first is not None:
        return set(range(1, min(first, total) + 1))
    elif last is not None:
        start = max(total - last + 1, 1)
        return set(range(start, total + 1))
    elif line_range is not None:
        try:
            a_str, b_str = line_range.split(",", 1)
            a, b = int(a_str), int(b_str)
            if a < 1 or b < a:
                raise ValueError
            return set(range(a, min(b, total) + 1))
        except ValueError:
            raise ValueError("Invalid -L format: expected 'a,b' with a<=b and a>=1")
    else:
        return set(range(1, total + 1))


def generate_arg_list(
    template: str,
    data: List[List[str]],
    first: int = None,
    last: int = None,
    line_range: str = None,
) -> Iterable[str]:
    """
    Generate output lines from Cartesian product of `data` using `template`.

    Parameters:
        template: string with placeholders like %1f%, %2f|1%, %i%
        data: list of lists; each inner list is a "column"
        first / last / line_range: same semantics as CLI options

    Yields:
        Expanded lines as strings.
    """
    if not data:
        return

    combos = list(product(*data))
    total = len(combos)
    if total == 0:
        return

    try:
        selected = select_indices(total, first, last, line_range)
    except ValueError as e:
        raise e

    for idx, combo in enumerate(combos, start=1):
        if idx in selected:
            yield substitute(template, list(combo), idx)


# def parse_args():
#     ap = argparse.ArgumentParser(add_help=False, prog="arg_list")
#     group = ap.add_mutually_exclusive_group()
#     group.add_argument(
#         "-f", type=int, metavar="N", help="output only first N combinations"
#     )
#     group.add_argument(
#         "-l", type=int, metavar="N", help="output only last N combinations"
#     )
#     group.add_argument(
#         "-L", metavar="a,b", help="output only combinations a‑b (inclusive, 1‑based)"
#     )
#     ap.add_argument("template", help="expansion template string")
#     ap.add_argument(
#         "-h", "--help", action="help", help="show this help message and exit"
#     )
#     return ap.parse_args()


def main():
    # Example in-memory data — replace as needed
    data: List[List[str]] = [
        generate_number_strings(1, 3),
        ["alice@example.com;Engineer", "bob@example.com;Designer"],
    ]

    try:
        # for line in generate_arg_list(
        #     template=args.template,
        #     data=data,
        #     first=args.f,
        #     last=args.l,
        #     line_range=args.L,
        # ):
        #     print(line)
        lines = generate_arg_list(
            template="%1f%: %2f;1%",
            data=data,
        )
        print("\n".join(lines))
    except ValueError as e:
        sys.stderr.write(f"arg_list: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
