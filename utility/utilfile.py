def read_list(filename):
    lst = []
    with open(filename) as file:
        for line in file:
            line = line.rstrip()
            if '_' not in line:
                if ',' in line:
                    line = line.split(',')[1].strip()
                try:
                    line = float(line)
                except ValueError:
                    pass
                lst.append(line)
    return lst


def write_list(lst, filename):
    file = open(filename, 'w')
    for value in lst:
        file.write(str(value) + '\n')
    file.close()
