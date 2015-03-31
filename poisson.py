#!usr/bin/python

class Node:
    def __init__(self, name, dist, nuc, left, right):
        self.name = name
        self.dist = dist
        self.nuc = nuc
        
        def setupArr(self):
            if self.nuc=='A':return [1,0,0,0]
            elif self.nuc == 'T':return [0,1,0,0]
            elif self.nuc == 'C':return [0,0,1,0]
            elif self.nuc == 'G':return [0,0,0,1]
            else: return [-1,-1,-1,-1]
        
        self.arr = setupArr(self)
        self.left = left
        self.right = right
    def __str__(self, depth=0):
            ret = ''
            tabs = '                     '*depth
            ret += tabs+'name: %s\n'%self.name+tabs+'branch length:%s\n'%self.dist+tabs+'nucleotide:%s\n'%self.nuc+tabs+'probs:%s'%self.arr+'\n'
            if self.left != None:
                ret += tabs +'Left:\n'+ self.left.__str__(depth+1)
            if self.right != None:
                ret += tabs+'Right:\n'+self.right.__str__(depth+1)
            return ret

    def likelihood(self):
        if self.left == None and self.right == None:
            pass
        else:
            for i in range(4):
                left = 0
                right = 0
                for j in range(4):
                    left += substitution_prob(i,j,self.left.dist)*self.left.arr[j]
                    right += substitution_prob(i,j,self.right.dist)*self.right.arr[j]
                self.arr[i] = left*right






tree2 = Node('5',None,'X',
                Node('4',15,'X',
                     Node('1',5,'A',None,None),
                     Node('2',10,'T',None,None)
                     ),
                Node('3',20,'C',None,None)
            )

def preorder(node):
    if node is not None:
        preorder(node.left)
        preorder(node.right)
        node.likelihood()
        print tree2

def substitution_prob(i,j,t):
    if i == j: #no net substitution
        return 1
    elif i+j == 3: #transition (assuming A=0, T=1, C=2, G=3
        return -t
    else: #transversion
        return t

preorder(tree2)
