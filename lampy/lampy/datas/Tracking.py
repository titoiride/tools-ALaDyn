from ..compiled_cython.read_tracking import total_tracking_read
import matplotlib.pylab as plt
import numpy as np
from ..utilities.Utility import speed_of_light, _tracking_directory, \
    _tracking_dictionary, _convert_component_to_index, _nearest_particle
import os

m_e = 0.511
e = 1.6E-7

class Tracking(object):

    def __new__(cls, Simulation):

        if not Simulation._tracking:
            print("Tracking not available")
            return None
        else:
            return super(Tracking, cls).__new__(cls)

    def __init__(self, Simulation):

        self._Simulation = Simulation
        self._params = Simulation.params
        self._directories = Simulation.directories
        self._timesteps = Simulation._timesteps
        self._dimensions = Simulation._dimensions
        self._path = Simulation.path
        self.normalized = False
        self.comoving = False
        self._Directories = Simulation._Directories
        self._stored_tracking = dict()
        self.iter_dictionary = dict()
        self._available_iterations = list()
        self._available_times = list()
        self._available_species = list()
        self._tracking = self._Simulation._tracking
        self._read_iter_dic()
        self._tracking_dictionary = _tracking_dictionary
        if self._params['n_dimensions'] == 2:
            self._array_dimensions = 7
        elif self._params['n_dimensions'] == 3:
            self._array_dimensions = 9
    
    def _check_index_in_ps(self, phase_space, index):

        comp = _convert_component_to_index( self._params, 'index')
        if index == 'all':
            return True
        mask = np.isin(phase_space[comp], index, assume_unique=True)
        return mask

    def _read_iter_dic(self):

        if not self._tracking:
            return

        for key, value in _tracking_dictionary.items():
            track_dic_path = os.path.join(self._path, \
                _tracking_directory, value)

            if os.path.isfile(track_dic_path):
                self._available_species += [key]
                self.iter_dictionary[key] = dict()
            elif not os.path.isfile(track_dic_path):
                continue

            with open(track_dic_path) as fp:
                lines = fp.readlines()

            # Removing header
            lines.pop(0)
            for line in lines:
                self._available_iterations += [int(line.split()[0])]
                self._available_times += [float(line.split()[1])]
                self.iter_dictionary[key][float(line.split()[1])] =\
                    int(line.split()[0])

    
    def _return_tracking_phase_space(self, timestep, species):

        track_list = self._search_track_by_timestep(timestep, species)
        if track_list is not None:
            pass
        else:
            self._track_read(timestep, species)


    def get_data(self, species, time):
        """
        Method that retrieves any field data and returns the 2D or 3D array.

        Parameters
        --------
        field_name : str
            Name of the array that is to be collected
        time : float
            Instant at which array is retrieved

        Returns
        --------
        field : dict
            Data are returned as a dictionary.
            field['data'] is the array containing the data
            field['time'] is the corresponding time
        """
        f = dict()
        time = self._Simulation._nearest_tracking_time(time, species)
        self._return_tracking_phase_space(time, species)
        f['data'] = self._stored_tracking[(time, species)]
        f['time'] = time

        return f

    def select_index(self, time=None, species=1, **kwargs):
        """
        Function that takes in input a phase space and 
        selects particles according to the given conditions.

        Parameters
        --------
        phase_space : dict or str
            If it is a string, is the phase space name.
            Otherwise, a phase space dictionary (i.e. a phase space collected
            via a get_data()) can be given.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        time : float, optional
            Time variable is needed when phase_space is a string

        Kwargs
        --------
        List of possible kwargs:

                'gamma_min', 'gamma_max', 'x_min', 'x_max', 'y_min', 'y_max',
                'z_min', 'z_max', 'weight_min', 'weight_max', 'nearest',
                'x', 'y', 'z'

                It is possible to select parts of phase space to analyze via
                the input kwargs.

        Returns
        --------
        ps_selected : np.array
            Phase space with the selected particles

        """
        n_dimensions = self._params['n_dimensions']

        possible_kwargs = ['gamma_min', 'gamma_max', 'x_min', 'x_max',
                           'y_min', 'y_max', 'z_min', 'z_max',
                           'weight_min', 'weight_max']

        var_dic_min = dict()
        var_dic_max = dict()
        var_dic_min['gamma_min'] = 'gamma'
        var_dic_max['gamma_max'] = 'gamma'
        var_dic_min['x_min'] = 'x'
        var_dic_min['y_min'] = 'y'
        var_dic_min['z_min'] = 'z'
        var_dic_max['x_max'] = 'x'
        var_dic_max['y_max'] = 'y'
        var_dic_max['z_max'] = 'z'
        var_dic_min['weight_min'] = 'weight'
        var_dic_max['weight_max'] = 'weight'
        
        if 'z_min' in kwargs and n_dimensions == 2:
            del kwargs['z_min']
        if 'z_max' in kwargs and n_dimensions == 2:
            del kwargs['z_max']

        nearest_part = False
        if 'nearest' in kwargs:
            if kwargs['nearest']:
                nearest_part = True
            if 'x' not in kwargs and 'y' not in kwargs and 'z' not in kwargs:
                print("""
        Error, when you ask 'nearest' you should specify the component
        'x', 'y' or 'z'
                    """)
                return
        comp_dict = dict()
        if 'x' in kwargs:
            comp = _convert_component_to_index(self._params, 'x')
            comp_dict[comp] = kwargs['x']
        if 'y' in kwargs:
            comp = _convert_component_to_index(self._params, 'y')
            comp_dict[comp] = kwargs['y']
        if 'z' in kwargs:
            comp = _convert_component_to_index(self._params, 'z')
            comp_dict[comp] = kwargs['z']
        comp_dict['index'] = _convert_component_to_index(self._params, 'index')

        index = list()
        kw_list = list(set(possible_kwargs).intersection(kwargs.keys()))

        time = self._Simulation._nearest_tracking_time(time, species)
        self._return_tracking_phase_space(time, species)

        ps = self._stored_tracking[(time, species)].copy()

        if nearest_part:
            ind = _nearest_particle(ps, comp_dict)
            return ind

        for kw in kw_list:
            if kw in var_dic_min.keys():
                comp = _convert_component_to_index(self._params, var_dic_min[kw])
                index.append(np.where(ps[comp] > kwargs[kw])[0])
            elif kw in var_dic_max.keys():
                comp = _convert_component_to_index(self._params, var_dic_max[kw])
                index.append(np.where(ps[comp] < kwargs[kw])[0])

        tot_index = index[0]
        for i in range(len(index)):
            tot_index = set(tot_index).intersection(index[i])
        tot_index = list(tot_index)

        index_comp = _convert_component_to_index(self._params, 'index')

        return np.array(ps[index_comp][tot_index])
        
    def scatterplot(self, time=None, component1='x', component2='y', species=1,
                index='all', comoving=False, s=1, **kwargs):
        """
        Method that produces a scatter plot of the given tracked phase space.

        It takes in input the tracked species and the time of interest.

        Parameters
        --------
        
        """
        if time is None:
            print("""
        Time not known, plotting last available dataset.
                """)
            time = self._available_times[-1]
        else:
            time = self._Simulation._nearest_tracking_time(time, species)

        self._return_tracking_phase_space(time, species)
        ps = self._stored_tracking[(time, species)].copy()

        mask = self._check_index_in_ps(ps, index)
        component1 = _convert_component_to_index(self._params, component1)
        component2 = _convert_component_to_index(self._params, component2)

        plt.scatter(ps[component1][mask], ps[component2][mask], s=s, **kwargs)

    def _search_track_by_timestep(self, timestep, species):

        track_list = None
        for key in self._stored_tracking.keys():
            if (timestep, species) == key:
                track_list = key

        return track_list

    def _track_read(self, timestep, species):

        from ..utilities.Utility import _sort_tracked_particles
        sim = self._Simulation
        ndim = self._params['n_dimensions']


        file_path = sim._derive_tracking_file_path(timestep, species)
        try:
            ps, part_number = total_tracking_read(file_path, self._params)
        except FileNotFoundError:
            print("""
        Tracking at time {} not available, impossible to read.
            """.format(timestep))
            raise

        index_in = _convert_component_to_index(self._params, 'index')

        ps_sorted = _sort_tracked_particles(ps, index_in)

        self._stored_tracking[(timestep, species)] = ps_sorted

    def trajectory_plot(self, time=None, component1='x', component2='y',
            species=1, index='all', **kwargs):

        from matplotlib.collections import LineCollection

        time = self._Simulation._nearest_tracking_time(time, species)
        time_index = self._available_times.index(time)
        time_array = self._available_times[:time_index + 1]
        max_time = max(time_array)
        min_time = max(time_array)
        plot_list_abscissa = dict()
        plot_list_ordinate = dict()

        if component1 != 'time':
            component1 = _convert_component_to_index(self._params, component1)
        component2 = _convert_component_to_index(self._params, component2)
        ind_comp = _convert_component_to_index(self._params, 'index')

        cmap_name = 'copper'
        norm = plt.Normalize(min_time, max_time)
        # Performing the first iteration explicitly to generate the numpy
        # arrays to be plotted 
        instant = time_array.pop(0)
        self._return_tracking_phase_space(instant, species)
        self._return_tracking_phase_space(instant, species)
        ps = self._stored_tracking[(instant, species)].copy()
        mask = self._check_index_in_ps(ps, index)
        if component1 != 'time':
            temp_comp1 = ps[component1][mask]
        temp_comp2 = ps[component2][mask]
        temp_index = ps[ind_comp][mask]
        particle_number = np.count_nonzero(mask)
        if component1 == 'time':
            for i in range(particle_number):
                plot_list_abscissa[temp_index[i]] = np.array([instant])
                plot_list_ordinate[temp_index[i]] = np.array(temp_comp2[i])
        else:
            for i in range(particle_number):
                plot_list_abscissa[temp_index[i]] = np.array(temp_comp1[i])
                plot_list_ordinate[temp_index[i]] = np.array(temp_comp2[i])

        for instant in time_array:
            self._return_tracking_phase_space(instant, species)
            ps = self._stored_tracking[(instant, species)].copy()
            mask = self._check_index_in_ps(ps, index)
            if component1 != 'time':
                temp_comp1 = ps[component1][mask]
            temp_comp2 = ps[component2][mask]
            temp_index = ps[ind_comp][mask]
            particle_number = np.count_nonzero(mask)
            for i in range(particle_number):
                if component1 == 'time':
                    plot_list_abscissa[temp_index[i]] = \
                        np.append(plot_list_abscissa[temp_index[i]], [instant])
                else:
                    plot_list_abscissa[temp_index[i]] = \
                        np.append(plot_list_abscissa[temp_index[i]], temp_comp1[i])
                plot_list_ordinate[temp_index[i]] = \
                    np.append(plot_list_ordinate[temp_index[i]], temp_comp2[i])

        for key in plot_list_abscissa.keys():
            max_xplot = max(plot_list_abscissa[key])
            max_yplot = max(plot_list_ordinate[key])
            min_xplot = min(plot_list_abscissa[key])
            min_yplot = min(plot_list_ordinate[key])

            # points = np.array([plot_list_abscissa[key],
            #                   plot_list_ordinate[key]])

            # segments = np.concatenate([points[:-1], points[1:]], axis=1)
            plt.plot( plot_list_abscissa[key], plot_list_ordinate[key])