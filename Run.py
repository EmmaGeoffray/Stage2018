import sys
from Definitions import *


############################################
# Run
############################################


# Maybe add help function later.

if os.path.isfile("./kappa_map.fits"):
    if os.path.isfile("./gamma_map.fits"):
        if not(os.path.isfile("./data_cube.fits")):
            # Tested default and non default
            # os.system("python Make_data_cube.py default")
            os.system("python Make_data_cube.py")
        else:
            print("Using pre-computed data_cube.fits.")

        test_path("./Maps")

        if len(sys.argv) > 1 and sys.argv[1] == "default":
            # Tested for default flux value
            # (have to check that there is no pbm with GERLUMPH for real flux values)
            cluster = Cluster()
            source = Source(mass=100*u.M_sun, z=2)
            n_pts = 15

            test_path("./Plots/M"+str(source.get_mass().value)+"_z"+str(source.get_z()))

            pix_for_adv_stats = cluster.basic_stats(source)
            cluster.adv_stats(pix_for_adv_stats, n_pts, source)

        elif len(sys.argv) > 1 and sys.argv[1] == "perso":
            # TO TEST
            z_near = 2
            z_med = 5
            z_far = 10
            z_farther = 20

            sources_near = [Source(1*u.M_sun, z_near), Source(10*u.M_sun, z_near), Source(30*u.M_sun, z_near)]
            sources_med = [Source(1*u.M_sun, z_med), Source(10*u.M_sun, z_med), Source(30*u.M_sun, z_med), Source(100*u.M_sun, z_med)]
            sources_far = [Source(1*u.M_sun, z_far), Source(10*u.M_sun, z_far), Source(30*u.M_sun, z_far), Source(100*u.M_sun, z_far), Source(300*u.M_sun, z_far)]
            sources_further = [Source(1*u.M_sun, z_farther), Source(10*u.M_sun, z_farther), Source(30*u.M_sun, z_farther), Source(100*u.M_sun, z_farther), Source(300*u.M_sun, z_farther)]

            cluster = Cluster()
            n_pts = 15

            print("Sources near (z={}):".format(z_near))
            for i in range(len(sources_near)):
                test_path("./Plots/M"+str(sources_near[i].get_mass().value)+"_z"+str(sources_near[i].get_z()))
                pix_for_adv_stats = cluster.basic_stats(sources_near[i])
                cluster.adv_stats(pix_for_adv_stats, n_pts, sources_near[i])

            print("Sources med (z={}):".format(z_med))
            for i in range(len(sources_med)):
                test_path("./Plots/M"+str(sources_med[i].get_mass().value)+"_z"+str(sources_med[i].get_z()))
                pix_for_adv_stats = cluster.basic_stats(sources_med[i])
                cluster.adv_stats(pix_for_adv_stats, n_pts, sources_med[i])

            print("Sources far (z={}):".format(z_far))
            for i in range(len(sources_far)):
                test_path("./Plots/M"+str(sources_far[i].get_mass().value)+"_z"+str(sources_far[i].get_z()))
                pix_for_adv_stats = cluster.basic_stats(sources_far[i])
                cluster.adv_stats(pix_for_adv_stats, n_pts, sources_far[i])

            print("Sources further (z={}):".format(z_further))
            for i in range(len(sources_further)):
                test_path("./Plots/M"+str(sources_further[i].get_mass().value)+"_z"+str(sources_further[i].get_z()))
                pix_for_adv_stats = cluster.basic_stats(sources_further[i])
                cluster.adv_stats(pix_for_adv_stats, n_pts, sources_further[i])

    else:
        print("Error: No gamma map.")
else:
    print("Error: No kappa map.")
