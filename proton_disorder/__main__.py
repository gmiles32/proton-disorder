from src import *

if __name__ == "__main__":

    parser = parser()
    args = parser.parse_args()
    dipole_target = args.dipole_target

    # Parse file, generate first configuration, test if it's good
    ice, box_dim = parse_input(args.input_file)
    h, hinv = gen_matrices(box_dim)
    ice = neighbour_list(ice,h,hinv)
    ice = init_hydrogens(ice)
    ice = adjust_bonds(ice)
    dipole = get_dipole(ice,h,hinv)

    # Optimize structure
    nround = 0
    while True:
        ice_old = ice
        dipole_old = dipole

        ice = shake_bonds(ice)
        ice = adjust_bonds(ice)

        dipole = get_dipole(ice,h,hinv)
        if (dipole < dipole_old):
            print("Round {}, dipole is now {:.5f}".format(nround,dipole))
        else:
            ice = ice_old
            dipole = dipole_old
        
        if (dipole < dipole_target):
            break
            
        nround += 1

    output(ice,h,hinv,filename=args.output_file,tip4p=(not args.tip3p))