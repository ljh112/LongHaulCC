import os
import argparse

def scale_cdf_file(input_file, output_file, scale_factor):
    """
    读取CDF分布表文件，将第一列的值扩大指定倍数
    
    参数:
    input_file -- 输入文件路径
    output_file -- 输出文件路径
    scale_factor -- 扩大倍数
    """
    # 读取原始文件
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # 处理每一行数据
    scaled_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2:
            try:
                value = float(parts[0]) * scale_factor
                percentage = parts[1]
                scaled_lines.append(f"{int(value) if value.is_integer() else value} {percentage}")
            except ValueError:
                # 如果无法转换为数字，保持原样
                scaled_lines.append(line.strip())
        else:
            # 不符合预期格式的行保持原样
            scaled_lines.append(line.strip())
    
    # 写入新文件
    with open(output_file, 'w') as f:
        f.write('\n'.join(scaled_lines))
    
    print(f"文件已处理完成。原始文件: {input_file}, 输出文件: {output_file}, 放大倍数: {scale_factor}")

def process_multiple_scales(input_file, output_prefix, scales):
    """
    处理多个放大倍数，为每个倍数生成对应的文件
    
    参数:
    input_file -- 输入文件路径
    output_prefix -- 输出文件名前缀
    scales -- 放大倍数列表
    """
    filename, ext = os.path.splitext(output_prefix)
    
    for scale in scales:
        output_file = f"{filename}_{scale}x{ext}"
        scale_cdf_file(input_file, output_file, scale)

if __name__ == "__main__":
    # 命令行参数解析
    parser = argparse.ArgumentParser(description='对CDF表格文件进行放大处理')
    parser.add_argument('--input', '-i', default="AliStorage2019.txt", help='输入文件名')
    parser.add_argument('--output', '-o', default="AliStorage2019_scaled.txt", help='输出文件名前缀')
    parser.add_argument('--scales', '-s', type=int, nargs='+', default=[1000], 
                        help='放大倍数列表，例如: 100 500 1000 2000')
    args = parser.parse_args()
    
    # 设置文件路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, args.input)
    output_prefix = os.path.join(current_dir, args.output)
    
    # 处理多个缩放比例
    process_multiple_scales(input_file, output_prefix, args.scales)
