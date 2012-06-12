'''
Created on May 25, 2012

@author: dstrauss
'''


class A(object):
    def __init__(self,inpt):
        self.a = inpt
        self.b = 20
        

class B(A):
    def __init__(self,inpt):
        self.c = 40
        self.d = inpt
        super(B,self).__init__(inpt)
        


Q = A(5)
print Q.a