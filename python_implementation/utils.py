

class toy_graph:
    def __init__(self,attr=[],adj=[]):
        self.adj=adj
        self.attr=attr
        self.graph_size=len(self.adj)


    def getsize(self):
        self.graph_size=len(self.adj)
        return self.graph_size