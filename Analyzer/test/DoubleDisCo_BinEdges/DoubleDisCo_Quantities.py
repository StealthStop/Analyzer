import numpy as np

# A simple wrapper class to a nested python dictionary
# and makes for adding getting things from said dictionary
# easier. All checking done here.
class Quantities:

    def __init__(self):
        self.quantities = {}
        self.finalEdges = ()

    def add(self, name, d1, d2, q):

        if d1 not in self.quantities:
            self.quantities[d1] = {}

        if d2 not in self.quantities[d1]:
            self.quantities[d1][d2] = {}

        self.quantities[d1][d2][name] = q

    def get(self, name, d1=None, d2=None):

        if name == "edges":
            if d1 != None and d2 != None:
                return self.finalEdges
            else:
                payload = []
                for d1, d2s in self.quantities.items():
                    for d2, _ in d2s.items():
                        payload.append((d1,d2))
                return payload

        if d1 == None and d2 == None:
            payload = []
            for d1, d2s in self.quantities.items():
                for d2, q in d2s.items():
                    payload.append(q[name])
                    
            return np.array(payload)

        if d1 != None and d1 not in self.quantities:
            raise LookupError("Cannot find the key \"%s\" in the dictionary !"%(d1))

        if d2 != None and d2 not in self.quantities[d1]:
            raise LookupError("Cannot find the key \"%s\" in the dictionary !"%(d2))

        return self.quantities[d1][d2][name]
        
    def getFinal(self, name):
        return self.get(name, self.finalEdges[0], self.finalEdges[1])
