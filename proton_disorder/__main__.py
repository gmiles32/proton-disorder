from proton_disorder import *
import sys

if __name__ == "__main__":

    parser = parser()
    args = parser.parse_args()
    # Check arguments
    if args.input_file == '':
        print("Please enter a valid input file")
        sys.exit()
    if args.output_file == '':
        print("Writing output to ice.xyz\n")
    if args.tip3p:
        print("Using TIP3P water in output\n")
    else:
        print("Using TIP4P water in output\n")
    if args.nshakes == -1:
        print("Will calculate optimal nshakes\n")
    else: 
        print("Nshakes is {}".format(args.nshakes))
   
    dipole_target = args.dipole_target
    print("Target dipole moment is: {:.3}\n".format(dipole_target))
    
    # Parse file, generate first configuration, test if it's good
    ice, box_dim = parse_input(args.input_file)
    h, hinv = gen_matrices(box_dim)
    ice = neighbour_list(ice,h,hinv)
    print("Created neighbour list.")
    ice = init_hydrogens(ice)
    print("Initialized hydrogens.")
    ice = adjust_bonds(ice)
    dipole = get_dipole(ice,h,hinv)
    print("Initial bond adjustment complete.")
    print("Starting dipole moment is {}".format(dipole))

    if args.nshakes == -1:
        nshakes = int(len(ice[OXY_COORD_INDEX]) / 10)
        print("Nshakes is {}".format(nshakes))
    else: 
        nshakes = args.nshakes
        print("Nshakes is {}".format(nshakes))
    
    print("Starting optimization. Output of the current lowest structure will be written to '{}.chk'\n".format(args.output_file))
    
    # Optimize structure
    nround = 0
    while True:
        ice_old = ice
        dipole_old = dipole

        ice = shake_bonds(ice,args.nshakes)
        ice = adjust_bonds(ice)

        dipole = get_dipole(ice,h,hinv)
        if (dipole < dipole_old):
            print("Round {}, dipole is now {:.5f}".format(nround,dipole))
            output(ice,h,hinv,filename=args.output_file,tip4p=(not args.tip3p),chk=True)        
        else:
            ice = ice_old
            dipole = dipole_old
        
        if (dipole < dipole_target):
            break
        
        nround += 1

    output(ice,h,hinv,filename=args.output_file,tip4p=(not args.tip3p))