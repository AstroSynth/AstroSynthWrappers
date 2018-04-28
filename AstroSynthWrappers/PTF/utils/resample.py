import numpy as np
from tqdm import tqdm


def general_sine(t, a, f, p): 
    return a * np.sin(2*np.pi*f*t+p)


def ephem_sum(t, functions):
    def f(x): 
        return general_sine(x, *functions[0])
    val = f(t)
    for func in functions[1:]:
        def f(x):
            return general_sine(x, *func)
        func_add = f(t)
        val += func_add
    return val


def mk_ephem_params(l, amp_range=[0, 0.02], freq_range=[0.0008333, 0.01667], phase_range=[0, 2*np.pi]):
    funcs = list()
    amps = np.random.uniform(amp_range[0], amp_range[1], size=(l,))
    freqs = np.random.uniform(freq_range[0], freq_range[1], size=(l,))
    phases = np.random.uniform(phase_range[0], phase_range[1], size=(l,))
    for amp, freq, phase in zip(amps, freqs, phases):
        funcs.append([amp, freq, phase])
    return funcs


def time_sample_PTF(n, data, pulsator, mag, l=3, noise_range=[0.02, 0.2], noise=True):
    if pulsator:
        ephem_params = mk_ephem_params(l)
    else:
        ephem_params = [[0, 0, 0]]
    t = np.array(data[n][0])
    noise_scatter = np.random.uniform(noise_range[0], noise_range[1])
    pure_ephem = ephem_sum(t, ephem_params)
    if noise:
        ephem = pure_ephem + np.random.normal(mag, noise_scatter, size=(len(t),))
    else:
        ephem = pure_ephem
    return t, ephem


def sto_sample_PTF(df, collection, pfrac, start=0, stop=None, noise=True):
    for freq, amp, (ID, n) in tqdm(df.xget_lc(start=start, stop=stop), total=stop-start):
        pulsator = bool(np.random.choice([0, 1], 1, p=[1-pfrac, pfrac])[0])
        t, csample = time_sample_PTF(n, df, pulsator, np.mean(collection.find_one({"_id":ID}, {"mag":1, "_id":0})['mag']), noise=noise)
        post = {"cSample":list(csample), "cSampleClass":pulsator}
        collection.update({"_id":ID}, {"$set":post}, upsert=False)


    

