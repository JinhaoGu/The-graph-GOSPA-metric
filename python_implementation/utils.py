

class toy_graph:
    def __init__(self,attr=None,adj=None):
        self.adj=[] if adj is None else adj
        self.attr=[] if attr is None else attr
        self.graph_size=len(self.adj)


    def getsize(self):
        self.graph_size=len(self.adj)
        return self.graph_size