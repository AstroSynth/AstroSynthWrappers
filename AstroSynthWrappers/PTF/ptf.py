import pandas as pd
import numpy as np

from scipy.signal import lombscargle
from scipy.optimize import curve_fit

import os
import subprocess
from itertools import tee

from pymongo import MongoClient
from bson.objectid import ObjectId

class PTFAstroSL:
    def __init__(self, collection, dbname='AstronomyData', name="PTFData", nk=100):
        with open('MongoDB.log', 'w') as f:
            self.mongodb = subprocess.Popen(["mongod"], stdout=f)
        self.name = name
        self.dbname = dbname
        self.collectionname = collection
        self.client, self.db, self.collection = self.__open_db_connection__()
        self.ordered_cursor = None
        self.size = self.collection.count()
        self.nk = nk

        self.__buffer__ = dict()
        self.__si__ = 0
        self.__fill_data_buffer__(0)
        self.dbft = False

    def use_db_ft(self, use):
        if 'Frequency' in self.collection.find_one() and use:
            self.dbft = True
        elif 'Frequency' not in self.collection.find_one():
            raise EnvironmentError('FT not cached on database')

    def __get_ordered_cursor__(self):
        return self.collection.aggregate([{"$sort" : {"size": -1}}], allowDiskUse=True)

    def __fill_data_buffer__(self, n):
        self.__si__ = n
        for i, lightcurve in enumerate(self.collection.find({"numerical_index": {"$lt":self.__si__+self.nk, "$gte":self.__si__}}, {"Frequency": 0, "Amplitude": 0})):
            self.__buffer__[i] = pd.DataFrame(data=lightcurve)

    def __fill_data_buffer_around__(self, n):
        self.__si__ = int(np.floor(n-0.5*self.nk))
        for i, lightcurve in enumerate(self.collection.find({"numerical_index": {"$lt":self.__si__+self.nk, "$gte":self.__si__}}, {"Frequency": 0, "Amplitude": 0})):
            self.__buffer__[i] = pd.DataFrame(data=lightcurve)

    def __lc_present__(self, n):
        if n < self.size:
            if self.__si__ <= n < self.__si__+self.nk:
                return True
            else:
                return False
        else:
            raise KeyError('key n={} not in range {}'.format(n, self.size))

    def __get_target_buffer__(self, n):
        if self.__lc_present__(n):
            return self.__buffer__[n-self.__si__]
        else:
            self.__fill_data_buffer_around__(n)
            return self.__buffer__[n-self.__si__]

    def get_data_from_db(self, n):
        if isinstance(n, int):
            data = self.collection.find_one({"numerical_index": n})
        elif isinstance(n, ObjectId):
            data = self.collection.find_one({"_id":n})
        return pd.DataFrame(data=data)

    def __open_db_connection__(self):
        client = MongoClient()
        db = client[self.dbname]
        collection = db[self.collectionname]
        return client, db, collection

    def __split__(self, frame, split_length=365, keyCol='obsHJD'):
        last_date = frame[keyCol].iloc[0]
        last_index = 0
        split_crit = list()
        for i, date in enumerate(frame[keyCol][1:]):
            print(date, split_length+last_date, date > split_length+last_date)
            if date > split_length+last_date:
                print('Gab found at {}'.format(i))
                stop = i+1
                split_crit.append([last_index, stop])
                last_index = stop
                last_date = date
        parts = [frame.iloc[x[0]:x[1]] for x in split_crit]
        parts.append(frame.iloc[split_crit[-1][1]:])
        return parts
    
    def normalize(self, arr):
        return [(y/np.mean(arr))-1 for y in arr]
    
    def get_ft(self, n=0, s=500, lock=False, nymult=1, num=1):
        if not self.dbft:
            return self.__generate_ft__(n=n, s=s, lock=lock, nymult=nymult, num=num)
        else:
            return self.__query_ft__(n=n)

    def __query_ft__(self, n=0):
        data = self.collection.find_one({"numerical_index": n})
        return data['Frequency'], data['Amplitude'], (data['_id'], data['numerical_index'])

    def __generate_ft__(self, n=0, s=500, lock=False, nymult=1, num=1):
        time, flux, meta = self.get_lc(n=n)

        if not len(time) <= 1:
            avg_sample_rate = (max(time)-min(time))/len(time)
            if avg_sample_rate != 0:
                ny = 1/(2*avg_sample_rate)
                res = 1/(max(time)-min(time))
            else:
                ny = 1/24
                res = 0.01
            if lock:
                us = s
            else:
                us = int(10*ny/(res))
            flux = self.normalize(flux)
            start_freq = 0.1*res
            end_freq = nymult*ny
            total_range = end_freq-start_freq
            fts = np.zeros((num, 2, us))

            for i in range(num):
                sub_start = ((i/num)*total_range) + start_freq
                sub_end = (((i+1)/num)*total_range) + start_freq
                f = np.linspace(sub_start, sub_end, us)
                pgram = lombscargle(np.array(time), np.array(flux), f, normalize=True)
                fts[i][0] = f
                fts[i][1] = pgram
            return fts[:, 0, :], fts[:, 1, :], meta

        else:
            f = np.linspace(0, 1/24, s)
            fts = np.zeros((num, 2, s))
            for i in range(num):
                fts[i][0] = f
                fts[i][1] = np.zeros(s)
            return fts[:, 0, :], fts[:, 1, :], meta

    def xget_orderd_lc(self, stop=None):
        if self.ordered_cursor is None:
            self.ordered_cursor = self.__get_ordered_cursor__()
        if stop is None:
            stop = self.size
        self.ordered_cursor, cur = tee(self.ordered_cursor)
        for i, target in enumerate(cur):
            if i >= stop:
                break
            yield self.get_lc(n=target['numerical_index'])

    def xget_orderd_ft(self, stop=None, s=500, lock=False, nymult=1, num=1):
        if self.ordered_cursor is None:
            self.ordered_cursor = self.__get_ordered_cursor__()
        self.ordered_cursor, cur = tee(self.ordered_cursor)
        if stop is None:
            stop = self.size
        for i, target in enumerate(cur):
            if i >= stop:
                break
            yield self.get_ft(n=target['numerical_index'], s=s, lock=lock, num=num)

    def xget_lc(self, stop=None, start=0):
        if stop is None:
            stop = self.size
        for i in range(start, stop):
            yield self.get_lc(n=i)

    def xget_ft(self, stop=None, start=0, s=500, lock=False, nymult=1, num=1):
        if stop is None:
            stop = self.size
        for i in range(start, stop):
            yield self.get_ft(n=i, s=s, lock=lock, num=num)

    def get_lc(self, n=0, se=0, full=True):
        data = self.__get_target_buffer__(n)
        if not full:
            data = self.__split__(data)
            data = data[se]
        return data.obsHJD.tolist(), data.mag.tolist(), (data._id[0], data.numerical_index[0])

    def get_object(self, n=0):
        return self.__get_target_buffer__(n)

    def xget_object(self, start=0, stop=None):
        if stop is None:
            stop = self.size
        for i in range(start, stop):
            yield self.get_object(n=i)

    def cache_ft(self, s=500, lock=False, nymult=1, num=1):
        for freq, amp, meta in self.xget_ft(s=s, lock=lock, nymult=nymult, num=num):
            ID = meta[0]
            post = {"Frequency":freq.tolist(), "Amplitude":amp.tolist()}
            self.collection.update({"_id":ID}, {"$set":post}, upsert=False)


    def __len__(self):
        return self.size

    def __repr__(self):
        out = list()
        out.append('PTF Data Wrapper')
        out.append('----------------')
        out.append('Size: {}'.format(self.size))
        out.append('Database: {}'.format(self.dbname))
        out.append('Collection: {}'.format(self.collectionname))
        return '\n'.join(out)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.get_lc(n=key)
        else:
            raise TypeError("PTF index must be type int")

    def __del__(self):
        print('Shutting Down Mongo Server')
        self.mongodb.terminate()


if __name__ == '__main__':
    PTFTest = PTFAstroSL('PTFData')




