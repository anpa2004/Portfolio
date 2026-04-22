import numpy as np
from attitude import Attitude
class ControlEvent:
    def __init__(self):
        print('initialized')

    def orbit_scheduler(self,t0,tf,dt=1,K=1/180,P=1/6):
        """
        Define a list of events to include
        """
        
        def mode2ref(m):
            if m == 1:
                return 'sun'
            elif m ==2:
                return 'gmo'
            else:
                return 'nadir'

        # Finding position lists for full orbits
        lmo = Attitude('lmo')
        gmo = Attitude('gmo')
        
        # Creating time vector
        N = int(round((tf-t0)/dt))
        t_vec = t0 + dt*np.arange(N+1)
        r_vecs_lmo = np.zeros((N+1,3),dtype=float)
        r_vecs_gmo = np.zeros((N+1,3),dtype=float)
        mode = np.zeros((N+1),dtype=int)

        # Simulating orbit and determining the state
        for i in range(N):
            r_lmo = lmo.orbit.r_vec(t_vec[i])
            r_gmo = gmo.orbit.r_vec(t_vec[i])

            r_vecs_lmo[i] = r_lmo
            r_vecs_gmo[i] = r_gmo

            if r_lmo[1] > 0:
                # Sun-pointing
                mode[i] = 1
            else:
                if np.abs(np.arccos(np.dot(r_lmo,r_gmo)/(np.linalg.norm(r_lmo)*np.linalg.norm(r_gmo)))) < np.deg2rad(35):
                    # GMO-pointing
                    mode[i] = 2
                else:
                    # Nadir-pointing
                    mode[i] = 3
        mode[-1] = mode[-2]
        
        # Creating schedules
        diff_vector = np.diff(mode)
        switches = np.nonzero(diff_vector)[0]
        schedules = [{'ref':mode2ref(mode[0]),'t_start':t0,'t_end':t_vec[switches[0]],'K':K,'P':P}]

        for i in range(len(switches) - 1):
            command = {'ref':mode2ref(mode[switches[i] + 1]),'t_start':t_vec[switches[i]+1],'t_end':t_vec[switches[i+1]],'K':K,'P':P}
            schedules.append(command)
        schedules.append({'ref':mode2ref(mode[switches[-1] + 1]),'t_start':t_vec[switches[-1]+1],'t_end':t_vec[-1],'K':K,'P':P})

        return {'schedules':schedules,'mode':mode,'switches':switches}