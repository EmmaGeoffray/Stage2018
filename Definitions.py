############################################
# Define classes and functions
# used in computations
############################################

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import G
from astropy.constants import c


# Tested
def random_pos(scale=100):
    x, y = np.random.choice(scale, 2)
    return x, y


# Tested
def random_vel(scale):
    # Write better random velocity later
    # Probably need to introduce a z dependence to be general
    vx = 0
    vy = 0
    while vx == 0 and vy == 0:
        vx, vy = np.random.choice(scale, 2)*np.random.choice([+1, -1], 2)
    return vx, vy


# Tested
def compute_sigma_cr(z_l, z_s, cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
    d_ls = cosmology.angular_diameter_distance_z1z2(z_l, z_s)
    d_l = cosmology.angular_diameter_distance(z_l)
    d_s = cosmology.angular_diameter_distance(z_s)

    sigma_cr = c**2*d_s/(d_ls*d_l*4*np.pi*G)
    sigma_cr = sigma_cr.to(u.M_sun/u.kpc**2)
    return sigma_cr


# Tested
def distances_ratio(z_l, z_s, cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
    d_ls = cosmology.angular_diameter_distance_z1z2(z_l, z_s)
    d_s = cosmology.angular_diameter_distance(z_s)
    ratio = d_ls/d_s
    return ratio.value


# Tested
def read_pos_and_m(pos_path):
    f = open(pos_path+"lens_pos.dat", "r")
    lines = f.readlines()
    x = []
    y = []
    m = []
    for i in lines:
        x.append(float(i.split()[0]))
        y.append(float(i.split()[1]))
        m.append(float(i.split()[2]))
    f.close()
    return x, y, m


# All tested without t units
class Simulation(object):

    def __init__(self, dt=1, n_maps=1, navg=50, res=512):
        self._dt = dt  # Should have units of time in the end ...
        self._n_maps = n_maps
        self._navg = navg
        self._res = res
        return

    def get_dt(self):
        return self._dt

    def get_n_maps(self):
        return self._n_maps

    def get_navg(self):
        return self._navg

    def get_res(self):
        return self._res


# All tested without t, pos or vel units
# TO CHECK
class Source(object):

    def __init__(self, mass=100*u.M_sun, z=10, ml_exp=3.5, ml_fact=1):
        self._mass = mass
        self._z = z
        self._ml_exp = ml_exp
        self._ml_fact = ml_fact
        return

    def get_mass(self):
        return self._mass

    def get_z(self):
        return self._z

    def get_ml_coeff(self):
        return self._ml_exp, self._ml_fact

    # Tested.
    # Still need to review random pos and vel.
    @staticmethod
    def trajectory(simulation=Simulation(), n_pts=15, lenses_moving=True):
        # Source parameters will be needed when introducing velocity dependence on z.
        # Meanwhile this method should remain static.

        dt = simulation.get_dt()
        n_maps = simulation.get_n_maps()

        x = [-1]
        y = [-1]
        while x[-1] > 99 or x[-1] < 0 or y[-1] > 99 or y[-1] < 0:
            del x[1:len(x)]
            del y[1:len(y)]
            x[0], y[0] = random_pos(100)
            # Source pos is in [0,99]x[0,99] (GERLUMPH boundaries)
            vx, vy = random_vel(5)

            if lenses_moving is True:
                for i in range(1, n_maps):
                    x.append(x[i-1] + vx*dt)
                    y.append(y[i-1] + vy*dt)
            else:
                for i in range(1, n_pts):  # Number n_pts chosen arbitrarily, see what is best.
                    x.append(x[i-1] + vx*dt)
                    y.append(y[i-1] + vy*dt)
        return x, y

    # TO CHECK
    def luminosity(self):
        lum = self._ml_fact*(self._mass/u.M_sun)**self._ml_exp * u.L_sun
        lum = lum.to(u.L_sun)
        return lum

    def app_mag_pre_mu(self, cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
        lum = self.luminosity()
        d = cosmology.distmod(self._z)
        abs_mag_sun = 4.83 * u.mag
        abs_mag = abs_mag_sun - (2.5*np.log10(lum.value) * u.mag)
        app_mag = d + abs_mag
        return app_mag

    def limit_mu(self, cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
        app_mag = self.app_mag_pre_mu(cosmology)
        app_mag_lim = 27 * u.mag  # General "average" limit for HFF images and is OK with the filters used here.
        mu_lim = 10**((app_mag-app_mag_lim).value/2.5)
        return mu_lim

    def app_mag_post_mu(self, mu=10**4, cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
        app_mag_post_mu = self.app_mag_pre_mu(cosmology) - 2.5*np.log10(mu)*u.mag
        return app_mag_post_mu


# Tested
def mu_above_limit(map_path, source=Source(), cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
    hdul_mu = fits.open(map_path)
    mu_data = np.absolute(hdul_mu[0].data)
    above_mu_lim = False
    if np.amax(mu_data) > source.limit_mu(cosmology):
        above_mu_lim = True
    return above_mu_lim


# Tested but need to review random vel
def move_lenses(simulation, path):
    # cluster param will be needed when putting vel dependence on z

    dt = simulation.get_dt()
    n_maps = simulation.get_n_maps()

    x, y, m = read_pos_and_m(path+"0/")

    vx = np.ones(len(x))
    vy = np.ones(len(y))
    for i in range(len(x)):
        vx[i], vy[i] = random_vel(5)

    for i in range(1, n_maps):
        x = x + vx*dt
        y = y + vy*dt

        f = open(path+str(i)+"/lens_pos.dat", 'w')
        for j in range(len(x)):
            f.write(str(x[j]) + "  " + str(y[j]) + "  " + str(m[j]) + "\n")
    return


# Tested
def new_map(kappa, gamma, kappa_st, simulation, path, begin=True):
    ks = (kappa-kappa_st)/kappa  # smooth matter fraction

    navg = simulation.get_navg()
    res = simulation.get_res()

    if begin is True:
        print("1st map")
        command = "srun -p gpu ./lenser_gpu -lens_pos r " + " -k " + str(kappa) + " -g " + str(gamma) + " -ks " + str(ks) + " -navg " + str(navg) + " -res " + str(res)
        print(command)
        os.system(command)
    else:
        command = "srun -p gpu ./lenser_gpu -lens_pos c -k " + str(kappa) + " -g " + str(gamma) + " -ks " + str(ks) + " -navg " + str(navg) + " -res " + str(res)
        print(command)
        os.system(command)

    # Make this nicer later ...
    os.system("mv ./lens_pos.dat " + path + "lens_pos.dat")
    os.system("mv ./map.bin " + path + "map.bin")
    os.system("mv ./mapmeta.dat " + path + "mapmeta.dat")

    os.system("python N2mu_bin2fits.py " + path)
    return


# Tested
def create_maps(kappa, gamma, kappa_st, simulation, path):
    if not os.path.exists(path):
        for i in range(simulation.get_n_maps()):
            os.makedirs(path+str(i)+"/")
    print("Files created.")

    new_map(kappa, gamma, kappa_st, simulation, path+"0/")
    print("new_map function worked at initialisation.")

    if simulation.get_n_maps() > 1:
        move_lenses(simulation, path)
        print("move_lenses function worked.")

        for i in range(1, simulation.get_n_maps()):
            os.system("mv " + path + str(i) + "/lens_pos.dat ./")
            new_map(kappa, gamma, kappa_st, simulation, path+str(i)+"/", False)
        print("new_map function worked at for all directories.")
    return


# Tested
class Cluster(object):

    def __init__(self, z=0.5440, cube_file='./data_cube.fits', ml_exp=1, ml_fact=1.5):
        self._z = z
        hdul = fits.open(cube_file)
        data = hdul[0].data
        hdr = hdul[0].header
        self._theta_x = np.deg2rad(hdr['CDELT1'])
        self._theta_y = np.deg2rad(hdr['CDELT2'])
        self._kappa = data[0]  # value for z_s=inf
        self._gamma = data[1]  # value for z_s=inf
        self._flux = data[2]  # units of erg/s/cm**2
        hdul.close()
        self._ml_exp = ml_exp
        self._ml_fact = ml_fact
        return

    def get_z(self):
        return self._z

    def get_kappa(self):
        return self._kappa

    def get_gamma(self):
        return self._gamma

    def get_flux(self):
        return self._flux

    def get_ml_coeff(self):
        return self._ml_exp, self._ml_fact

    def get_pix_dim(self):
        return self._theta_x, self._theta_y

    # TO CHECK (for 20 order of magnitude problem)
    def compute_sigma_st(self, cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
        d = cosmology.angular_diameter_distance(self._z)

        # Small angle approximation.
        pix_area = np.absolute(self._theta_x)*d * np.absolute(self._theta_y)*d  # TO CHANGE using right header keywords and distance d.
        pix_area = pix_area.to(u.cm**2)
        # Flux units are in erg/s/cm**2.
        luminosity = self._flux*pix_area.value * u.erg / u.s
        sigma_st = (luminosity/u.L_sun/self._ml_fact)**(1/self._ml_exp) / pix_area * u.M_sun
        sigma_st = sigma_st.to(u.M_sun/u.kpc**2)
        return sigma_st

    # TO CHECK (for 20 order of magnitude problem)
    def compute_kappa_star(self, source=Source(), cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
        sigma_cr = compute_sigma_cr(self._z, source.get_z(), cosmology)
        sigma_st = self.compute_sigma_st(cosmology)
        kappa_st = (sigma_st/sigma_cr.to(sigma_st.unit)).value
        ratio = distances_ratio(self._z, source.get_z(), cosmology)
        if np.any(kappa_st > self._kappa*ratio):
            print("Error: kappa_star > kappa.")
        return kappa_st

    # Tested
    # TO CHECK (for 20 order of magnitude problem)
    def basic_stats(self, source=Source(), simul=Simulation(), cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
        print("Basic statistics for a source of mass {} and redshift {}.".format(source.get_mass(), source.get_z()))

        ratio = distances_ratio(self._z, source.get_z(), cosmology)

        # Need to multiply kappa and gamma by this ratio to get the values for a specific redshift of the source
        # (and not z_s=inf as we defined for the Cluster instance variables)
        kappa = self._kappa*ratio
        gamma = self._gamma*ratio
        kappa_st = self.compute_kappa_star(source, cosmology)
        one = np.ones(kappa.shape)
        mu = np.absolute(1/((one-kappa)**2-gamma**2))

        mu_thr = 10**(int(np.log10(source.limit_mu(cosmology)))-1)  # exponent is (order_of_magnitude-1)

        n_above_thr = 0
        n_above_lim = 0

        nx, ny = kappa.shape

        pix_for_adv_stats = []
        pix_above_lim = np.zeros((nx, ny))

        for i in range(nx):
            for j in range(ny):
                if mu[i, j] > mu_thr:
                    n_above_thr = n_above_thr+1
                    path = "./Maps/k" + str(kappa[i, j]) + "_g" + str(gamma[i, j]) + "_kst" + str(kappa_st[i, j]) + "/"
                    create_maps(kappa[i, j], gamma[i, j], kappa_st[i, j], simul, path)
                    if mu_above_limit(path+"0/map.fits", source, cosmology):
                        n_above_lim = n_above_lim+1
                        pix_above_lim[i, j] = 1
                        pix_for_adv_stats.append((i, j))

        percent_above_thr = 100*float(n_above_thr)/(nx*ny)
        percent_above_lim = 100*float(n_above_lim)/(nx*ny)
        print('Number of pixels above \mu_thr=' + str(mu_thr) + ' (for z_s=' + str(source.get_z()) + '): ' + str(n_above_thr))
        print('Percentage of pixels above \mu_thr=' + str(mu_thr) + ' (for z_s=' + str(source.get_z()) + '): ' + str(percent_above_thr) + '%')
        print('Number of pixels above \mu_lim=' + str(source.limit_mu(cosmology)) + ' (for z_s=' + str(source.get_z()) + '): ' + str(n_above_lim))
        print('Percentage of pixels above \mu_lim=' + str(source.limit_mu(cosmology)) + ' (for z_s=' + str(source.get_z()) + '): ' + str(percent_above_lim) + '%')

        d = cosmology.angular_diameter_distance(self._z)

        # Trigonometry using small angle approximation
        dx = np.absolute(self._theta_x)*d
        dy = np.absolute(self._theta_y)*d

        x = np.linspace(0, dx*nx, num=nx, endpoint=False)
        y = np.linspace(0, dy*ny, num=ny, endpoint=False)
        x_mesh, y_mesh = np.meshgrid(x, y)

        plt.pcolormesh(x_mesh.value, y_mesh.value, pix_above_lim, cmap=plt.cm.get_cmap('Blues', 2))
        plt.colorbar(ticks=[0, 1])
        plt.clim(0, 1)
        plt.title('Pixels above $\mu$ limit')
        plt.xlabel('x [{}]'.format(x_mesh.unit))
        plt.ylabel('y [{}]'.format(y_mesh.unit))
        plt.savefig("./Plots/M"+str(source.get_mass().value)+"_z"+str(source.get_z())+"/above_mu_map.png")
        plt.clf()
        return pix_for_adv_stats

    # Tested
    # TO CHECK: (for 20 order of magnitude problem)
    # TO CHECK: what are best conditions for good light curves ?
    # TO CHANGE (once know about Gerlumph distances)
    def adv_stats(self, pix_for_adv_stats, n_pts=15, source=Source(), simul=Simulation(), cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
        print("Advanced statistics.")

        ratio = distances_ratio(self._z, source.get_z())

        kappa = self._kappa*ratio
        gamma = self._gamma*ratio
        kappa_st = self.compute_kappa_star(source, cosmology)
        good_light_curves = np.zeros(self._kappa.shape)
        light_curve_example = False

        for n in range(len(pix_for_adv_stats)):
            i, j = pix_for_adv_stats[n]

            dir_map = "./Maps/k" + str(kappa[i, j]) + "_g" + str(gamma[i, j]) + "_kst" + str(kappa_st[i, j]) + "/0/map.fits"
            hdul_mu = fits.open(dir_map)
            mu_data = np.absolute(hdul_mu[0].data)

            for m in range(100):
                xs, ys = source.trajectory(simul, n_pts, lenses_moving=False)
                mu = np.ones(len(xs))
                mag = np.ones(len(xs))

                for i in range(n_pts):
                    mu[i] = mu_data[xs[i], ys[i]]
                    mag[i] = source.app_mag_post_mu(mu[i], cosmology).value

                # TO CHECK: what are best conditions for good light curve statistics ?
                mu_3rd_percentile = np.percentile(mu, 75)
                mag_diff = np.amax(mag)-np.amin(mag)

                if mu_3rd_percentile > source.limit_mu(cosmology) and mag_diff > 1:
                    good_light_curves[i, j] = good_light_curves[i, j]+1

                    if light_curve_example is False:
                        light_curve_example = True
                        angle = np.sqrt(np.power(xs-xs[0]*np.ones(len(xs)), 2)+np.power(ys-ys[0]*np.ones(len(ys)), 2))
                        ds = cosmology.angular_diameter_distance(source.get_z())

                        # TO CHANGE !!! Once know about Gerlumph distances
                        gerlumph_coeff = 1

                        # Trigonometry using small angle approximation
                        distance = ds * angle * gerlumph_coeff

                        plt.plot(distance, mag)
                        plt.title('Good light curve example')
                        # TO CHANGE !!! Once know about Gerlumph distances
                        # plt.xlabel('Distance [{}]'.format(distance.unit))
                        plt.xlabel('Distance [?]')
                        plt.ylabel('Magnitude [mag]')
                        plt.savefig("./Plots/M"+str(source.get_mass().value)+"_z"+str(source.get_z())+"/light_curve_example.png")
                        plt.clf()

        nx, ny = kappa.shape
        d = cosmology.angular_diameter_distance(self._z)

        # Trigonometry using small angle approximation
        dx = np.absolute(self._theta_x)*d
        dy = np.absolute(self._theta_y)*d

        x = np.linspace(0, dx*nx, num=nx, endpoint=False)
        y = np.linspace(0, dy*ny, num=ny, endpoint=False)
        x_mesh, y_mesh = np.meshgrid(x, y)

        plt.pcolormesh(x_mesh.value, y_mesh.value, good_light_curves, cmap=plt.cm.get_cmap('Blues', 5))
        plt.colorbar(ticks=[0, 20, 40, 60, 80, 100])
        plt.clim(0, 100)
        plt.title('Good light curves percentage')
        plt.xlabel('x [{}]'.format(x_mesh.unit))
        plt.ylabel('y [{}]'.format(y_mesh.unit))
        plt.savefig("./Plots/M"+str(source.get_mass().value)+"_z"+str(source.get_z())+"/good_light_curves_percentage.png")
        plt.clf()
        return


# Tested without t, pos or vel units
# TO CHANGE when have Gerlumph units
def plot_moving_lenses(kappa, gamma, kappa_st, source=Source(), simulation=Simulation(), cosmology=FlatLambdaCDM(H0=70, Om0=0.3)):
    dir_map = "./Maps/k" + str(kappa) + "_g" + str(gamma) + "_kst" + str(kappa_st) + "/"
    print("Moving source.")
    x, y = source.trajectory(simulation)
    mu = np.ones(len(x))
    mag = np.ones(len(x))

    print("Preparing plot:")
    for i in range(len(x)):
        print(i)
        hdul_mu = fits.open(dir_map+str(i)+"/map.fits")
        mu_data = np.absolute(hdul_mu[0].data)

        plt.imshow(mu_data, norm=LogNorm(), vmin=10**0, vmax=10**2)
        plt.colorbar(ticks=[10**0, 10**1, 10**2])
        # Think again about how to do this colorbar thing right
        plt.scatter(x[i], y[i], color='k')
        plt.title('$\mu$ map')
        plt.xlabel('x [?]')
        plt.ylabel('y [?]')
        plt.savefig(dir_map+"/map"+str(i)+".png")
        plt.clf()

        mu[i] = mu_data[x[i], y[i]]
        mag[i] = (source.app_mag_post_mu(mu[i], cosmology)).value

    d = np.sqrt(np.power(x-x[0]*np.ones(len(x)), 2)+np.power(y-y[0]*np.ones(len(y)), 2))

    plt.plot(d, mu)
    plt.title('Source magnification')
    plt.xlabel('Distance [?]')
    plt.ylabel('$\mu$ []')
    plt.savefig(dir_map+"/mu_source.png")
    plt.clf()

    plt.plot(d, mag)
    plt.title('Source magnitude')
    plt.xlabel('Distance [?]')
    plt.ylabel('Magnitude [mag]')
    plt.savefig(dir_map+"/mag_source.png")
    plt.clf()
    return d, mu, mag


# Tested
def test_path(path_folder):
    if os.path.isdir(path_folder):
        os.system("rm -rf " + path_folder + "/*")
        print("Removed all files in " + path_folder + " folder.")
    else:
        os.system("mkdir " + path_folder)
    return


############################################
# Tests
############################################

# simul = Simulation(dt=2,n_maps=3)
# print(simul.get_n_maps())
# print(simul.get_dt())
# print(simul.get_navg())
# print(simul.get_res())
# print(simul.get_ml_coeff())

# source = Source(mass=100*u.M_sun)
# simulation = Simulation(n_maps=1)
# print("The source is at redshift {} and weights {}.".format(source.get_z(), source.get_mass()))
# print(source.trajectory(simulation))

# source2 = Source()
# mu = source2.limit_mu()
# print("The limit magnification for a star of mass {} at redshift {} is {}.".format(source2.get_mass(), source2.get_z(), np.format_float_scientific(mu)))

# cluster = Cluster()
# print(cluster.get_z())
# print(cluster.get_kappa()[337, 332])
# print(cluster.get_gamma()[337, 332])
# print(cluster.get_flux()[337, 332])
# print(cluster.get_header())
# print(cluster.compute_sigma_st()[337, 332])

# cluster = Cluster()
# source = Source()
# print(cluster.compute_kappa_star(source)[337,332])

# simul = Simulation(n_maps=6)
# source = Source()
# kappa = 1.3
# gamma = 0.5
# kappa_st = 0.2
# path = "./Maps/k" + str(kappa) + "_g" + str(gamma) + "_kst" + str(kappa_st) + "/"
# create_maps(kappa, gamma, kappa_st, simul, path)
# plot_moving_lenses(kappa, gamma, kappa_st, source, simul)
