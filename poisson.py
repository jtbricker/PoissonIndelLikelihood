#!usr/bin/python

import copy, math
import numpy as np
from scipy.linalg import expm

class color: #Just for easy read formatting
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

class Node:
    def __init__(self, name, brLen, left=None, right=None):
        self.name = name
        self.brLen = brLen
        self.left = left
        self.right = right
        self.base = 'X'
        self.fhats = [0.0]*len(eps_alphabet)    #same order as eps_alphabet
        self.fhat = 0.0
        self.isLeaf = left == None and right == None
        self.p_surv = self.calc_psurv()  #Survival Probability given insertion at node
        self.f = 0.0
        self.ins_prior = 0.0
    
    def __str__(self, depth=0):
        ret = ''
        tabs = '                     '*depth
        ret += tabs + color.UNDERLINE + 'name: %s\n'%self.name + color.END
        ret += tabs + 'branch length:%s\n'%self.brLen
        ret += tabs + 'nucleotide:%s\n'%self.base
        ret += tabs + 'f_hats:%s\n'%self.fhats
        ret += tabs + 'f_hat:%s\n'%self.fhat
        ret += tabs + 'p_surv:%s\n'%self.p_surv
        ret += tabs + 'f: %s\n'%self.f
        ret += tabs + 'insertion prior: %s\n'%self.ins_prior
        ret += tabs + 'Leaf:%s'%self.isLeaf+'\n'
        if not self.isLeaf:
            ret += tabs + color.PURPLE + 'Left:\n'+ color.END + self.left.__str__(depth+1)
            ret += tabs+ color.GREEN + 'Right:\n' + color.END +self.right.__str__(depth+1)
        return ret

    def set_base(self, tname, tbase):
        if self.name == tname:
            self.base = tbase
            self.set_fhats()
        else:
            if not self.isLeaf:
                self.left.set_base(tname,tbase)
                self.right.set_base(tname,tbase)

    def set_fhats(self):
        if self.base == 'a':
            self.fhats = [1.0,0.0]
        elif self.base == '-':
            self.fhats = [0.0,1.0]
        elif self.base == 'X':
            self.fhats = [0.0,0.0]

    def calc_fhats(self):
        if self.isLeaf:
            temp_fhat=0.0
            for k in range(len(alphabet)):
                temp_fhat = temp_fhat+ self.fhats[k]*alpha_freq[k]
            self.fhat = temp_fhat
            return   #No need to calculate fhats at leaves
        for i in range(len(eps_alphabet)):  #To calculate each of fhats (for each base and gap)
            left_sum = 0   #inner sum for left child
            right_sum = 0   #inner sum for right child
            for j in range(len(eps_alphabet)):  #Loop over possible values of each child
                left_sum = left_sum+ self.subs_prob(j,i,self.left.brLen)*self.left.fhats[j]
                right_sum = right_sum+ self.subs_prob(j,i,self.right.brLen)*self.right.fhats[j]
            self.fhats[i] = left_sum*right_sum
        temp_fhat = 0.0
        for k in range(len(alphabet)):
            temp_fhat = temp_fhat+ self.fhats[k]*alpha_freq[k]
        self.fhat = temp_fhat
                    
    def subs_prob(self, from_base, to_base, t):
        S = expm(t*np.matrix(theta))
        return S[from_base,to_base]

    def calc_psurv(self):
        if self.brLen == 0:
            return 1 #For convenience
        else:
            return (1-math.exp(-mu*self.brLen))/(mu*self.brLen)

    def calc_f(self, num_extant):
        if self.brLen == 0:
            self.f = self.fhat
        else:
            if num_extant == 0: #Special calculation of f for null column
                self.f = 1 + self.p_surv*(self.fhat-1)
            else:
                if calc_num_extant(self) == num_extant:
                    self.f = self.p_surv * self.fhat
                else:
                    self.f = 0

    def calc_ins_prior(self):
        if self.brLen == 0:
            mfact = 1/mu
        else:
            mfact = self.brLen
        self.ins_prior = mfact/(tau + 1/mu)

class Column:
    def __init__(self, tree_template, align_col,site_num):
        self.tree = self.fill_tree(tree_template,align_col)
        self.site_num = site_num
        self.num_extant = calc_num_extant(self.tree)
        self.col_prob = 0.0
    
    def fill_tree(self, tree_template, align_col):
        for node in range(len(align_col[0])):
            tree_template.set_base(align_col[0][node],align_col[1][node])
        return tree_template
    
    def __str__(self):
        ret = color.RED + 'Site %s\n'%self.site_num + color.END
        ret += 'Number of Extant Leafs: %s\n'%self.num_extant
        ret += 'Column Probability: %s\n'%self.col_prob
        ret += '\n'
        ret += self.tree.__str__()
        return ret

def calc_num_extant(node):
    if node.isLeaf:
        if node.base != 'X' and node.base != '-':
            return 1
        else:
            return 0
    else:
        return calc_num_extant(node.left) + calc_num_extant(node.right)

def preorder_fhats(node):
    if node is not None:
        preorder_fhats(node.left)
        preorder_fhats(node.right)
        node.calc_fhats()

def site_stats(node, num_extant):    #calculate f values and insertion priors for all nodes
    if node.isLeaf:
        node.calc_f(num_extant)
        node.calc_ins_prior()
    else:
        node.calc_f(num_extant)
        node.calc_ins_prior()
        site_stats(node.left,num_extant)
        site_stats(node.right,num_extant)

def calc_total_brLen(node):
    if node.isLeaf:
        return node.brLen
    else:
        return node.brLen + calc_total_brLen(node.left) + calc_total_brLen(node.right)

def calc_col_prob(node):
    if node.isLeaf:
        return node.ins_prior*node.f
    else:
        return node.ins_prior*node.f + calc_col_prob(node.left) + calc_col_prob(node.right)

def phi(null_col_prob, num_col, v):
    return math.pow(v,num_col)*math.exp((null_col_prob-1)*v)/math.factorial(num_col)

def ptm(phi, columns): #marginal likelihood p_t(m)
    product = phi
    for column in columns:
        product = product*column.col_prob
    return product

if __name__ == '__main__':
    #GIVEN: Tree, Alignment, PIP Parameters, Substitution matrix, Alphabet, and


    lamb = 2.0 #insertion rate
    mu = 1.0 #deletion rate
    theta = [[-lamb ,mu],[lamb,-mu]] #substitution matrix
    alphabet = ['a']
    alpha_freq = [1]

    eps_alphabet = alphabet + ['-']  #Epsilon-augmented alphabet

    tree = Node('v1',0,
                 Node('v0',1,
                      Node('v2',1,None,None),
                      Node('v3',1,None,None)
                      ),
                 Node('v4',2,None,None)
                 )
    tau = calc_total_brLen(tree) #normalizaton of tree (sum of all tree branch lengths)
    v = lamb*(tau + (1/mu))
    
    alignment = [['v2','v3','v4'],['-','a','a'],['a','a','-']]
    null_align = ['-']*len(alignment[0])

    sites = [0]*(len(alignment)-1)
    for k in range(len(alignment)-1):
        sites[k] = Column(copy.deepcopy(tree), [alignment[0],alignment[k+1]],k+1)
    null_site = Column(copy.deepcopy(tree), [alignment[0],null_align],'Null')

    for site in (sites + [null_site]):
        preorder_fhats(site.tree)

    sites[0].tree.fhat = 0.012
    sites[1].tree.fhat = 0.043
    sites[1].tree.left.fhat = 0.14
    null_site.tree.fhat = 0.67

    for site in (sites + [null_site]):
        site_stats(site.tree, site.num_extant)
        site.col_prob = calc_col_prob(site.tree)
        print site

    vphi = phi(null_site.col_prob, len(sites), v)
    vptm = ptm(vphi, sites)
    print '\n\n\nMarginal Likelihood:  %s\n' %vptm
    print 'Log Marginal Likelihood:  %s\n' %math.log(vptm)







