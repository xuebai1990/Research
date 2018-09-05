fout = open("tip4p2005.data", 'w')

nmol = 7180
boxlength = 60.0

fout.write("Molecule TIP4P/2005 water\n")
fout.write("\n")
fout.write(str(nmol*3)+" atoms\n")
fout.write(str(nmol*2)+" bonds\n")
fout.write(str(nmol)+" angles\n")
fout.write("\n")
fout.write("2 atom types\n")
fout.write("1 bond types\n")
fout.write("1 angle types\n")
fout.write("\n")
fout.write("0.0 "+str(boxlength)+"      xlo xhi\n")
fout.write("0.0 "+str(boxlength)+"      ylo yhi\n")
fout.write("0.0 "+str(boxlength)+"      zlo zhi\n")
fout.write("\n")
fout.write("Atoms\n")
fout.write("\n")
counter = 0
atom_counter = 0
mol_counter = 0
with open("config1a.dat", 'r') as f:
    for line in f:
        counter += 1
        if counter % 2 == 0 or counter % 8 == 7:
            continue
        else:
            atom_counter += 1
            if atom_counter % 3 == 1:
                atom_type = 2
            else:
                atom_type = 1
            mol_counter = (atom_counter - 1) // 3 + 1
            numbers = line.split()
            if atom_counter % 3 == 1:
                numbers.insert(0, "-1.1128")
            else:
                numbers.insert(0, "0.5564")
            numbers.insert(0, str(atom_type))
            numbers.insert(0, str(mol_counter))
            numbers.insert(0, str(atom_counter))
            numbers.append('\n')
            fout.write(' '.join(numbers))

    fout.write("\n")
    fout.write("Bonds\n")
    fout.write("\n")
    bond_counter = 0
    for i in range(mol_counter):
        bond_counter += 1
        numbers = [str(bond_counter), '1']
        numbers.append(str(3*i+1))
        numbers.append(str(3*i+2))
        fout.write(' '.join(numbers)+'\n')
        bond_counter += 1
        numbers = [str(bond_counter), '1']
        numbers.append(str(3*i+1))
        numbers.append(str(3*i+3))
        fout.write(' '.join(numbers)+'\n')

    fout.write("\n")
    fout.write("Angles\n")
    fout.write("\n")
    angle_counter = 0
    for i in range(mol_counter):
        angle_counter += 1
        numbers = [str(angle_counter), '1']
        numbers.extend([str(3*i+2), str(3*i+1), str(3*i+3)])
        numbers.append("\n")
        fout.write(' '.join(numbers))



