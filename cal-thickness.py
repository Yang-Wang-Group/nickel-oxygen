import glob

O_TYPE = 3    
NI_TYPE = 6   
OUTPUT_FILE = 'analysis_results.csv'  
Z_THRESHOLD = 75 

def safe_float_format(value, default='NaN'):
    
    try:
        return f"{float(value):.4f}"
    except (ValueError, TypeError):
        return default

def parse_dump_file(filename):
    with open(filename, 'r') as f:
       
        oxygen_z_above = []  
        oxygen_z_below = [] 
        nickel_z_above = []  
        nickel_z_below = []  
        
        
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith('ITEM: ATOMS'):
                headers = line.strip().split()[2:]
                try:
                    type_idx = headers.index('type')
                    z_idx = headers.index('z')
                except ValueError:
                    print(f"警告：{filename} 缺少type或z列")
                    break
                
               
                for line in f:
                    if line.startswith('ITEM:'): 
                        break
                    parts = line.strip().split()
                    if len(parts) <= max(type_idx, z_idx):
                        continue
                    try:
                        atom_type = int(parts[type_idx])
                        z = float(parts[z_idx])
                        
                        if atom_type == O_TYPE:
                            if z > Z_THRESHOLD:
                                oxygen_z_above.append(z)
                            else:
                                oxygen_z_below.append(z)
                        elif atom_type == NI_TYPE:
                            if z > Z_THRESHOLD:
                                nickel_z_above.append(z)
                            else:
                                nickel_z_below.append(z)
                    except (ValueError, IndexError):
                        continue
                break
    
    return oxygen_z_above, oxygen_z_below, nickel_z_above, nickel_z_below


file_list = glob.glob('dump.pos.*')
file_list.sort(key=lambda f: int(f.split('.')[-1]))

with open(OUTPUT_FILE, 'w') as out_f:
    
    out_f.write('Filename,O_min_z_above_75,Ni_max_z_above_75,O_max_z_below_75,Ni_min_z_below_75\n')
    
    for filename in file_list:
        try:
            o_above, o_below, ni_above, ni_below = parse_dump_file(filename)
            
            
            o_min_above = min(o_above) if o_above else None
            ni_max_above = max(ni_above) if ni_above else None
            o_max_below = max(o_below) if o_below else None
            ni_min_below = min(ni_below) if ni_below else None
            
            
            o_min_above_str = safe_float_format(o_min_above)
            ni_max_above_str = safe_float_format(ni_max_above)
            o_max_below_str = safe_float_format(o_max_below)
            ni_min_below_str = safe_float_format(ni_min_below)
            
            
            out_f.write(f"{filename},{o_min_above_str},{ni_max_above_str},{o_max_below_str},{ni_min_below_str}\n")
        except Exception as e:
            print(f"处理文件 {filename} 时出错: {str(e)}")
            out_f.write(f"{filename},NaN,NaN,NaN,NaN\n")

print(f"处理完成！结果已保存到 {OUTPUT_FILE}")