# fct_formatter.py
import sys

def IP_to_id(ip):
    return int((int(ip,16)-0x0b000001)/0x00000100)

def format_fct_line(line):
    parts = line.strip().split()
    if len(parts) != 8:
        return f"# Invalid line: {' '.join(parts)}"
    
    # 将时间字段从ns转为ms（1ms=1e6ns）
    ns_to_ms = lambda x: round(float(x)/1e6, 6)  # 保留三位小数

    fields = {
        "Source Node": IP_to_id(parts[0]),
        "Destination Node": IP_to_id(parts[1]),
        "Source Port": parts[2],
        "Destination Port": parts[3],
        "Traffic Size(Bytes)": parts[4],
        "Start Time(ms)": ns_to_ms(parts[5]),
        "Multiple Flows Completion Time(ms)": ns_to_ms(parts[6]),
        "Single Flow Completion Time(ms)": ns_to_ms(parts[7])
    }
    
    return "\n".join([f"{k}: {v}" for k, v in fields.items()])

def main():
    if len(sys.argv) != 2:
        print("Usage: python fct_formatter.py fct_ls.txt")
        return
    
    filename = sys.argv[1]
    with open(filename, 'r') as f:
        for i, line in enumerate(f, 1):
            formatted = format_fct_line(line)
            print(f"----- Line {i} -----\n{formatted}\n")

if __name__ == "__main__":
    main()
# cd /home/sdn/ljh/LongHaulCC/simulation/mix
# python fct_formatter.py /home/sdn/ljh/LongHaulCC/simulation/mix/fct_ls.txt