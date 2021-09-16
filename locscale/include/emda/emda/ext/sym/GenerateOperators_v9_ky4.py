#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:51:10 2021

@author: garib
"""
import numpy as np
import random as rnd
import math as math
import sys
eps_l = 1.0e-5


def from_two_axes_to_group_v2(axis1,order1,axis2=None, order2=0,toler=1.0e-2):
    #
    # This routine  acccepts two axes and corresponding two orders.
    #
    # Here we assume that we have at least one of the major axes.
    # I.e. for any group the highest order axis is given. 
    # For dihedral - n≥2
    # For T - n = 3, for O n=4 and for I n=5
    #
    # Generating from other combinations these axes should be done in anohter routine.
    #
    if axis2 is None: axis2 = np.array([0.0,0.0,1.0])
    axes_out = list(); grp_out = list()
    order_out = list()
    grp_symbol = ' '
    if order2 <= 0 :
        grp_symbol = 'C'+str(order1)
        grp_now = generate_cyclic(axis1, order1)
        axes_out.append(axis1); order_out.append(order1); grp_out = grp_now
        return(grp_symbol,order_out,axes_out,grp_out)
    maxorder = max(order1,order2)
    minorder = min(order1,order2)
    cangle = np.dot(axis1,axis2)/(np.linalg.norm(axis1)*np.linalg.norm(axis2))
    #
    # is min order == 2 and angle between axis 90 (pi/2)
    if minorder == 2 and np.abs(cangle)< toler :
        grp_symbol = 'D'+str(maxorder)
        #
        #  Make sure that axes are orthogonal. a la Gram-Schmidt orthogonolisation. Just in case. 
        if order1 == maxorder :
            axis_l1 = axis1; order_l1 = order1
            axis_l2 = axis2; order_l2 = order2
        else :
            axis_l1 = axis2; order_l2 = order2
            axis_l2 = axis1; order_l2 = order1
        axis_l1 = axis_l1/np.linalg.norm(axis_l1)
        axis_l2 = axis_l2 - np.dot(axis_l1,axis_l2)/np.dot(axis_l1,axis_l1)*axis_l1
        axis_l2 = axis_l2/np.linalg.norm(axis_l2)
        order_out,axes_out,grp_out = generate_all_elements(axis_l1,order_l1,axis_l2,order_l2,toler)
        return(grp_symbol,order_out,axes_out,grp_out)
    #
    #  Remaining groups are T, O and I
    #
    # Is it T?
    if maxorder == 3 :
        if minorder == 2 and (np.abs(cangle - 0.5773) < toler or np.abs(cangle+0.5773)<toler) :
            grp_symbol = 'T'
        elif minorder == 3 and (np.abs(cangle+0.3333) < toler or np.abs(cangle-0.3333) < toler) :
            grp_symbol = 'T'
        #order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        #if grp_symbol != ' ' :
        #    return(grp_symbol,order_out,axes_out,grp_out)
    elif maxorder == 4 :
        if minorder == 4 and np.abs(cangle) < toler :
            grp_symbol = 'O'
        if minorder == 3 and (np.abs(cangle+0.5773) or np.abs(cangle+0.5773)<toler):
            grp_symbol='O'
        #order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        #if grp_symbol != ' ' :
        #    return(grp_symbol,order_out,axes_out,grp_out)  
    elif maxorder == 5 :
        if minorder == 5 and (np.abs(cangle+0.4472)< toler or np.abs(cangle-0.4472)<toler):
            grp_symbol = 'I'
            #order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        if minorder == 3 and (np.abs(cangle+0.7947)< toler or np.abs(cangle-0.7947) < toler 
                              or np.abs(cangle-0.187592)<toler or np.abs(cangle+0.187592)<toler):
            grp_symbol='I'
            #order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        if minorder == 2 and (np.abs(cangle-0.85065) < toler or np.abs(cangle+0.85065) < toler or 
                              np.abs(cangle-0.52373) < toler or np.abs(cangle+0.52373)< toler):
            grp_symbol='I'
            #order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        #if(grp_symbol != ' '):
        #    return(grp_symbol,order_out,axes_out,grp_out)
    if grp_symbol != ' ':
        order_out,axes_out,grp_out = generate_all_polyhedron(grp_symbol,axis1,order1,axis2,order2)
        return(grp_symbol,order_out,axes_out,grp_out)
    #
    # IF we have not returned yet then raise an error
    print('A finite group could not be identified')
    grp_symbol = ' '
    sys.exit('A finite group could not be identified')
    #return(grp_symbol,order_out,axes_out,grp_out)

def generate_all_polyhedron(grp_symbol,axis1,order1,axis2,order2):
    if grp_symbol == 'T':
        order_out,axes_out,grp_out = generate_all_tetrahedron(axis1,order1,axis2,order2)
    elif grp_symbol == 'O':
        order_out,axes_out,grp_out = generate_all_octahedron(axis1,order1,axis2,order2)
    elif grp_symbol == 'I':
        order_out,axes_out,grp_out = generate_all_icosahedron(axis1,order1,axis2,order2)
    else :
        print('Group symbol must be one of: T, O, I')
        sys.exit('Wrong group symbol: must be one of: T,O, I')
    return(order_out,axes_out,grp_out)
    pass

def from_two_axes_to_group(axis1,order1,axis2 = np.array([0.0,0.0,1.0]),order2=0,toler=1.0e-4):
    #
    # This routine  acccepts two axes and corresponding two orders.
    #
    # Here we assume that we have at least one of the major axes.
    # I.e. for any group the highest order axis is given. 
    # For dihedral - n≥2
    # For T - n = 3, for O n=4 and for I n=5
    #
    # Generating from other combinations these axes should be done in anohter routine.
    #
    axes_out = list(); grp_out = list()
    order_out = list()
    grp_symbol = ' '
    if order2 <= 0 :
        grp_symbol = 'C'+str(order1)
        grp_now = generate_cyclic(axis1, order1)
        axes_out.append(axis1); order_out.append(order1); grp_out = grp_now
        return(grp_symbol,order_out,axes_out,grp_out)
    maxorder = max(order1,order2)
    minorder = min(order1,order2)
    cangle = np.dot(axis1,axis2)/(np.linalg.norm(axis1)*np.linalg.norm(axis2))
    #
    # is min order == 2 and angle between axis 90 (pi/2)
    if minorder == 2 and np.abs(cangle)< toler :
        grp_symbol = 'D'+str(maxorder)
        order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2,toler)
        return(grp_symbol,order_out,axes_out,grp_out)
    #
    #  Remaining groups are T, O and I
    #
    # Is it T?
    if maxorder == 3 :
        if minorder == 2 and (np.abs(cangle - 0.5773) < toler or np.abs(cangle+0.5773)<toler) :
            grp_symbol = 'T'
        elif minorder == 3 and (np.abs(cangle+0.3333) < toler or np.abs(cangle-0.3333) < toler) :
            grp_symbol = 'T'
        order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        if grp_symbol != ' ' :
            return(grp_symbol,order_out,axes_out,grp_out)

    if maxorder == 4 :
        if minorder == 4 and np.abs(cangle) < toler :
            grp_symbol = 'O'
        if minorder == 3 and (np.abs(cangle+0.5773) or np.abs(cangle+0.5773)<toler):
            grp_symbol='O'
        order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        if grp_symbol != ' ' :
            return(grp_symbol,order_out,axes_out,grp_out)
    
    if maxorder == 5 :
        if minorder == 5 and (np.abs(cangle+0.4472)< toler or np.abs(cangle-0.15224)<toler 
                              or np.abs(cangle-0.4472)<toler):
            grp_symbol = 'I'
            order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        if minorder == 3 and (np.abs(cangle+0.7947)< toler or np.abs(cangle-0.7947)<toler):
            grp_symbol='I'
            order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        if minorder == 2 and (np.abs(cangle-0.8507) < toler or np.abs(cangle+0.8507) < toler):
            grp_symbol='I'
            order_out,axes_out,grp_out = generate_all_elements(axis1,order1,axis2,order2)
        if(grp_symbol != ' '):
            return(grp_symbol,order_out,axes_out,grp_out)
    #
    # IF we have not returned yet then raise an error
    print('A finite group could not be identified')
    grp_symbol = ' '
    sys.exit('A finite group could not be identified')
    #return(grp_symbol,order_out,axes_out,grp_out)
    
#import line_profiler
#profile = line_profiler.LineProfiler()
#import atexit
#atexit.register(profile.print_stats)
#@profile
def generate_all_elements(axis1,order1,axis2 = np.array([0.0,0.0,1.0]),order2=0,toler=1.0e-6) :
    # 
    # Generate all group elements. Output will be as a list of cyclic groups.
    #
    grp_out = list()
    order_out = list()
    axes_out = list()
    if order2 == 0 :
        grp_out = generate_cyclic(axis1,order1)
        #grp_out.append(grp_n)
        order_out.append(order1)
        axes_out.append(axis1)
        return(order_out,axes_out,grp_out)
    grp_out = generate_cyclic(axis1,order1)
    order_out.append(order1); axes_out.append(axis1)
    grp1 = generate_cyclic(axis2,order2)
    grp_out = add_groups_together(grp_out,grp1)
    order_out.append(order2); axes_out.append(axis2)
    
    grp_out_new = copy_group(grp_out)
    things_to_do = True
    while things_to_do :
        things_to_do = False
        for i, ri in enumerate(grp_out) :
            for j, rj in enumerate(grp_out) :
                if i !=0 and j!=0 and i != j:
                    r3 = np.dot(ri,rj)
                    if not is_in_the_list_rotation(r3,grp_out) :
                        things_to_do = True
                        order3,axis3 = find_order(r3)
                        grp3 = generate_cyclic(axis3,order3)
                        grp_out_new = add_groups_together(grp_out_new,grp3)
                        order_out.append(order3)
                        axes_out.append(axis3)
                    r3 = np.dot(rj,ri)
                    if not is_in_the_list_rotation(r3,grp_out) :
                        things_to_do = True
                        order3,axis3 = find_order(r3)
                        grp3 = generate_cyclic(axis3,order3)
                        grp_out_new = add_groups_together(grp_out_new,grp3)
                        order_out.append(order3)
                        axes_out.append(axis3)
        grp_out = copy_group(grp_out_new)
    #
    # Filter out axes (if they are parallel to each other then select one, for the order we should take
    # the highest order). In our case it should happen only for the group O. We may have 2 and four fold symmetries
    # with the same axis
    axes_out_new = list(); order_out_new = list()
    for i,axisi in enumerate(axes_out):
        order_cp = order_out[i]
        for j,axisj in enumerate(axes_out) :
            if i < j:
                cangle = np.dot(axisi,axisj)/(np.linalg.norm(axisi)*np.linalg.norm(axisj))
                if np.abs(cangle-1.0) < toler or np.abs(cangle+1) < toler :
                    # same axis
                    if order_cp <= order_out[j] :
                        order_cp= 0
                        break
        if order_cp > 0 :
            axes_out_new.append(axisi); order_out_new.append(order_cp)
    return(order_out_new,axes_out_new,grp_out)

def copy_group(grp):
    grp_out = list()
    for i,r in enumerate(grp):
        grp_out.append(r)
    return(grp_out)
    
def find_order(r,toler=1.0e-3) :
    order=1
    r_id = np.identity(3)
    r3 = np.copy(r_id)
    things_to_do = True
    A = r_id
    while things_to_do and order < 100:
        things_to_do = False
        r3 = np.dot(r3,r)
        if np.sum(np.abs(r3-r_id)) > toler :
            A = A + r3
            things_to_do = True
            order = order + 1
        else :
            continue
    if order >= 100 :
        sys.exit("The order of the group is too high: order > 100")
    A = A/order
    axis_l = find_axis(A)
        
    return(order,axis_l)
            
            
def add_groups_together(grp_in,grp_add):
    grp_out = list()
    for i,r in enumerate(grp_in):
        grp_out.append(r)
    #grp_out = np.copy(grp_in)
    for i,r in enumerate(grp_add) :
        if not is_in_the_list_rotation(r,grp_out) :
            grp_out.append(r)
    return(grp_out)
    
def generate_cyclic(axis,order):
    #
    #  This function generates all cyclic group elements using axis and order of the group
    if order <=0 or np.sum(np.abs(axis)) < eps_l :
        print("order = ",order, "axis = ",axis)
        sys.exit("Either order or axis is zero")
    gout = list()
    id_matr = np.identity(3)
    gout.append(id_matr)
    angle = 2.0*np.pi/order
    axis = axis/np.linalg.norm(axis)
    exp_matr = np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    axis_outer = np.outer(axis,axis)
    m_int = id_matr - axis_outer
    for i in range(order-1):
        angle_l = angle*(i+1)
        stheta = np.sin(angle_l); ctheta = np.cos(angle_l)
        m_l = exp_matr*stheta + m_int*ctheta +axis_outer
        gout.append(m_l)
    return(gout)
        

def AngleAxis2rotatin(axis,angle):
    #
    #  Convert axis and ange to a rotation matrix. Here we use a mtrix form of the relatiionship
    # IT may not be the moost efficient algorithm, but it should work (it is more elegant)
    if np.sum(np.abs(axis)) < eps_l :
        print("Axis = ",axis," Angle = ",angle)
        sys.exit("Axis is zero")
    id_matr = np.identity(3)
    axis = axis/np.sqrt(np.dot(axis,axis))
    exp_matr = np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    axis_outer = np.outer(axis,axis)
    m_int = id_matr - axis_outer
    stheta = np.sin(angle); ctheta = np.cos(angle)
    m_l = exp_matr*stheta + m_int*ctheta +axis_outer
    return(m_l)

def Rotation2AxisAngle_cyclic(m_in,eps_l = 1.0e-5):
    #
    # Here we assume that rotation matrix is an element of a cyclic group
    # This routine gives the smallest angle for this cyclic group. 
    # To find axis of the rotation we use the fact that if we define 
    # A = 1/n sum_i-0^(n-1) (R^i) then this operator is a projector to the axis of rotation
    # i.e. for Ax will be on hte axis for any x. IT could be equal 0, in this case we select another x
    A = m_in
    m1 = m_in
    id_matr = np.identity(3)
    cycle_number = 1
    ended = False
    while not ended and cycle_number < 200:
        if np.sum(np.abs(m1-id_matr)) < eps_l :
            ended = True
            break
        m1 = np.dot(m1,m_in)
        A = A + m1
        cycle_number = cycle_number + 1
    # take a ranom vector
    if cycle_number >= 150 :
        print("matrix ",m_in)
        print("Try to change the tolerance: eps_l = XXX")
        sys.exit("The matrix does not seem to be producing a finite cyclic group")
    A = A/cycle_number
    axis = np.array([0.0,0.0,0.0])
    while np.dot(axis,axis) < eps_l:
        xin = np.array([rnd.random(),rnd.random(),rnd.random()])
        axis = np.dot(A,xin)
        if np.dot(axis,axis) > eps_l :
            axis = axis/np.sqrt(np.dot(axis,axis))
    if axis[2] < 0.0 :
        axis = -axis
    elif axis[2] == 0.0 and axis[1] < 0.0:
        axis = -axis
    angle = 2.0*np.pi/cycle_number
    return(axis,angle,cycle_number)
    
def Rotation2AxisAngle_general(m_in):
    #
    #  This routine shuld work for any rotation matrix
    axis = np.array([1,0.0,0.0])
    angle = math.acos(max(-1.0,np.min((np.trace(m_in)-1)/2.0)))
    if np.sum(np.abs(m_in-np.transpose(m_in))) < eps_l :
        # 
        # It is a symmetric matrix. so I and m_in form a cyclic group
        A = (np.identity(3) + m_in)/2.0
        axis = np.array([0.0,0.0,0.0])
        while np.linalg.norm(axis) < eps_l :
            axis = np.dot(A,np.array([rnd.random(),rnd.random(),rnd.random()]))
    else :
        axis[0] = m_in[1,2] - m_in[2,1]
        axis[1] = m_in[0,2] - m_in[2,0]
        axis[2] = m_in[0,1] - m_in[1,0]
    if axis[2] < 0.0 :
        axis = -axis
        angle = 2.0*np.pi - angle
    elif axis[2] < eps_l and axis[1] < 0.0 :
        axis = -axis
        angle = 2.0*np.pi - angle
    axis = axis/np.linalg.norm(axis)
    return(axis,angle)

def is_it_cyclic(m_in,eps_l = 1.0e-4):
    #
    # Find out if this matrix is a member of a finite cyclic group. 
    # If yes return order of the group and the axis
    id_matr = np.identity(3)
    found = False
    m1 = m_in
    A = m1
    cycle_number = 1
    while not found and cycle_number < 100:
        if np.sum(np.abs(m1-id_matr)) < eps_l:
            found = True
            break
        m1 = np.dot(m1,m_in)
        A = A + m1
        cycle_number = cycle_number + 1
    A = A/cycle_number
    axis = np.array([0,0.0,0.0])
    while np.linalg.norm(axis) < eps_l :
        axis = np.dot(A,np.array([rnd.random(),rnd.random(),rnd.random()]))
    if axis[2] < 0.0 :
        axis = -axis
    elif axis[2] == 0.0 and axis[1] < 0.0 :
        axis = -axis
    axis = axis/np.linalg.norm(axis)
    return(found,axis,cycle_number)

#@profile
def is_in_the_list_rotation(m_in,m_list,toler = 1.0e-3):
    id_matr = np.identity(3)
    return np.any(np.abs(np.trace(np.dot(np.transpose(m_in), m_list)-id_matr[:,None], axis1=0,axis2=2)) < toler)

def closest_rotation(m_in,m_list):
    id_matr = np.identity(3)
    return min(np.abs(np.trace(np.dot(np.transpose(m_in), m_list)-id_matr[:,None], axis1=0,axis2=2)))

def is_it_same_axis(axis1,axis2,eps_l = 1.0e-5) :
    output = False
    if np.dot(axis1,axis2)/(np.linalg.norm(axis1)*np.linalg.norm(axis2)) > 1-eps_l:
        output = True
    return(output)

#
#  Simple routines to find greatest common divisors (gcd) and least common multiples (lcm)
def find_gcd(n1,n2):
    x,y = n1,n2
    if n1 < n2 :
        x,y=n2,n2
    if x==y:
        return(x)
    
    while (y):
        x,y = y,x%y
    return(x)

def find_lcm(n1,n2):
    g = find_gcd(n1,n2)
    return(n1*n2/g)

def two_axes_to_third(order1,order2,info=False):
    theta1 = 2.0*np.pi/order1; theta2 = 2.0*np.pi/order2
    ct1 = np.cos(theta1); st1 = np.sin(theta1)
    ct2 = np.cos(theta2); st2 = np.sin(theta2)
    a1 = ct1*ct2-ct1-ct2+1; b1 = st1*st2; c1 = ct1*ct2+ct1+ct2
    r1 = np.identity(3)
    r1[0,0] = ct1; r1[0,1]=-st1; r1[1,0] = st1; r1[1,1]=ct1
    #
    if info :
        print('---------------------------------------------------------------------------------')
        print('----  We must remember that only certain combination of axes is possible     ----')
        print('----  Cyclic groups may have any combination of orders of rotation as sooon  ----')
        print('----  axes coincide                                                          ----')
        print('----  Dihedral groups are generated by n-fold and two-fold symmetries if     ----')
        print('----  their axes are orthogonal each other.                                  ----')
        print('----  Remaining possible groups are:                                         ----')
        print('----  I - icosahedral group: it has 2,3 and 5-fold axes                      ----')
        print('----  O - octohedral group : it has 2,3 and 4-fold axes                      ----')
        print('----  T - tetraheral group : it has 2 and 3-fold axes                        ----')
        print('----  No other combination of symmetry axes is possible                      ----')
        print('---------------------------------------------------------------------------------')
    for i in range(4):
        order3 = i+2
        theta3 = 2.0*np.pi/order3
        ct3 = np.cos(theta3)  #; st3 = np.sin(theta3)
        c_l = c1 - 2.0*ct3-1
        D = b1**2 - a1*c_l
        if D < 0.0:
            print(' oder1 ', order1, ' and order2 ',order2,' does not generate any order ',order3,' rotation axis')
        elif abs(D) < 1.0e-8:
            z = -b1/a1
            z = np.round(z,8)
            #if z < 0.0 : 
            #    z = -z
            if( z >= -1.0 and z <=1.0) :
                angle = math.acos(max(-1,min(1.0,z)))*180.0/np.pi
                print('cos of the angle between first two axes = ',z)
                print_axes_combination(order1,order2,order3,angle)
                y=np.sqrt(1.0-z**2)
                r2 = AngleAxis2rotatin(np.array([0.0,y,z]), theta2)
                r3 = np.dot(r1,r2)
                #
                #  Calculate the projector on the axis
                axis_1 = find_axis_projector(r3,order3)
                print('Axis 1: ',np.array([0.0,0.0,1.0]))
                print('Axis 2: ',np.array([0.0,y,z]))
                print('Axis 3: ',axis_1)
        else :
            z1 = (-b1+np.sqrt(D))/a1; z2 = (-b1-np.sqrt(D))/a1
            if z1 >= -1.0 and z1 <= 1.0 :
                angle = math.acos(z1)*180.0/np.pi
                print('cos of the angle between first two axes = ',z1)
                print_axes_combination(order1,order2,order3,angle)
                y = np.sqrt(min(1.0,1-z1**2))
                r2 = AngleAxis2rotatin(np.array([0.0,y,z1]), theta2)
                r3 = np.dot(r1,r2)
                axis_1 = find_axis_projector(r3,order3)
                print('Axis 1: ',np.array([0.0,0.0,1.0]))
                print('Axis 2: ',np.array([0.0,y,z1]))
                print('Axis 3: ',axis_1)
            if z2 >= -1.0 and z2 <= 1.0 :
                angle = math.acos(z2)*180.0/np.pi
                print('cos of the angle between first two axes = ',z2)
                print_axes_combination(order1,order2,order3,angle)
                y = np.sqrt(min(1.0,1-z2**2))
                r2 = AngleAxis2rotatin(np.array([0.0,y,z2]), theta2)
                r3 = np.dot(r1,r2)
                axis_1 = find_axis_projector(r3,order3)
                print('Axis 1: ',np.array([0.0,0.0,1.0]))
                print('Axis 2: ',np.array([0.0,y,z2]))
                print('Axis 3: ',axis_1)
                
def find_axis(amatr):
    
    #
    #  We assume that amatr is a projector. I.e. y = amatr x is on the the symmetry axis. 
    # To avoid problem of 0 vector we try several times to make sure that 0 vector is not generated
    axis1 = np.array([0.0,0.0,0.0])
    while np.linalg.norm(axis1) < 0.001 :
        axis1 = np.dot(amatr,np.array([rnd.random(),rnd.random(),rnd.random()]))
    axis1 = axis1/np.linalg.norm(axis1)
    axis1 = np.around(axis1,10)
    axis1 = axis1/np.linalg.norm(axis1)
    #axis1 = axis1[0]
    # Remove annoying negative signs
    #np.where(np.abs(axis1)<1.0e-8,0.0,axis1)
    axis1[axis1==0.]=0.
    
    if axis1[2] < 0.0 :
        axis1 = -axis1
    elif axis1[2] == 0.0 and axis1[1] < 0.0:
        axis1 = -axis1
    return(axis1)

def find_axis_projector(r,order):
    
    #
    #  We assume that amatr is a projector. I.e. y = amatr x is on the the symmetry axis. 
    # To avoid problem of 0 vector we try several times to make sure that 0 vector is not generated
    amatr = calculate_projector_on_the_axis(r,order)
    axis1 = np.array([0.0,0.0,0.0])
    while np.linalg.norm(axis1) < 0.001 :
        axis1 = np.dot(amatr,np.array([rnd.random(),rnd.random(),rnd.random()]))
    axis1 = axis1/np.linalg.norm(axis1)
    axis1 = np.around(axis1,10)
    axis1 = axis1/np.linalg.norm(axis1)
    #
    # Remove annoying negative signs
    #np.where(np.abs(axis1)<1.0e-8,0.0,axis1)
    axis1[axis1==0.]=0.
    return(axis1)
def  print_axes_combination(order1,order2,order3,angle) :
    print(order1,'-fold and ', order2, '-fold axes can generate ',order3,
          '-fold axis if the angle between them is', angle)
    
def calculate_projector_on_the_axis(r,order):
    amatr = np.identity(3)
    r4 = np.identity(3)
    for i in range(order-1):
        r4 = np.dot(r4,r)
        amatr = amatr + r4
    amatr = amatr/order
    return(amatr)
        
def find_all_axes(matrices,toler=1.0e-4):
    #
    # For a given set of matrices this function returns axes and corresponding orders
    # It is assumed that each matrix is a member of a finite cyclic group
    axes = list()
    orders = list()
    for i,r in enumerate(matrices):
        order1,axis1 = find_order(r)
        axes.append(axis1)
        orders.append(order1)
    
    axes_out_new = list(); order_out_new = list()
    for i,axisi in enumerate(axes):
        order_cp = orders[i]
        for j,axisj in enumerate(axes) :
            if i < j:
                cangle = np.dot(axisi,axisj)/(np.linalg.norm(axisi)*np.linalg.norm(axisj))
                if np.abs(cangle-1.0) < toler or np.abs(cangle+1) < toler :
                    # same axis
                    if order_cp <= orders[j] :
                        order_cp= 0
                        break
        if order_cp > 0 :
            axes_out_new.append(axisi); order_out_new.append(order_cp)
    return(order_out_new,axes_out_new)
    
def find_minimal_axis(axes_in,order_in) :
    if len(order_in) == 1:
        axis1,order1 = axes_in[0],order_in[0]
        return(axis1,order1,np.array([0,0,1]),0)
    axis2,order2 = np.array([0,0,0]),1
    return(axis1,order1,axis2,order2)

def operators_from_symbol(op):
    op_l = op.upper()
    if op_l == 'I' :
        order1 = 5
        axis1 = np.array([0, 0.85065, 0.52573])
        order2 = 5
        axis2 = np.array([0.52573, 0, 0.85065])
    elif op_l == 'O' :
        order1 = 4
        axis1 = np.array([0,0,1.0])
        order2 = 4
        axis2 = np.array([0,1.0,0])
    elif op_l[0] == 'T' :
        order1 = 3
        axis1 = np.array([0.0,0.0,1.0])
        order2 = 3
        axis2 = np.array([0.0,0.94280904,-0.33333333])
    elif op_l[0] == 'D' :
        order1 = int(op[1:])
        axis1 = np.array([0,0,1.0])
        order2 = 2
        axis2 = np.array([1.,0.,0.])
    elif op_l[0] == 'C' :
        order1 = int(op[1:])
        axis1 = np.array([0,0,1.0])
        order2 = 0
        axis2 = np.array([0,0,0.0])
    return(generate_all_elements(axis1,order1,axis2,order2))

def generate_all_tetrahedron(axis1,order1,axis2,order2,toler=1.0e-1):
    orders,axes = tetrahedron_axes()
    axis_l1 = axis1/np.linalg.norm(axis1); axis_l2 = axis2/np.linalg.norm(axis2)
    if order1 == 3 and order2 == 3:
        found = 0
        for i,o in enumerate(orders):
            if o==3 and found !=2 :
                if "ax1" in locals() :
                    ax2 = axes[i]
                    found = found + 1
                else:
                    ax1 = axes[i]
                    found = found + 1
    elif order1 == 3 and order2 == 3:
        #
        # Any two three fold axes would do
        found = 0
        for i,o in enumerate(orders):
            if (o == 3 or o == 2)  and found != 2:
                if o == 3:
                    ax1 = axes[i]
                    found = found + 1
                else:
                    ax2 = axes[i]
                    found = found + 1
    else:
        #
        # Error
        print('Only combination of two and three or three and three folds are allowed')
        sys.exit('Wrong combinations of axes orders: must be (3,3) or (2,2)')
    ax1 = ax1/np.linalg.norm(ax1)
    ax2 = ax2/np.linalg.norm(ax2)
    #
    #  Find transformation matrix. Make sure that angles between axes is correct. Make sure that we
    #  consider ax2 and -ax2, just in case. We assume that the lengths of vectors are equal to 1.
    cangle = np.dot(axis_l1,axis_l2)
    cangle_t = np.dot(ax1,ax2)
    if np.abs(cangle+cangle_t) < toler :
        axis_l2 = -axis_l2
    ax3 = np.cross(ax1,ax2); ax3 = ax3/np.linalg.norm(ax3)
    axis_l3 = np.cross(axis_l1,axis_l2); axis_l3 = axis_l3/np.linalg.norm(axis_l3)
    AA = np.array([axis_l1,axis_l2,axis_l3])
    BB = np.array([ax1,ax2,ax3])
    from scipy.linalg import orthogonal_procrustes
    R,sca = orthogonal_procrustes(AA,BB)
    print("Rg: ", R)
    axes_out = list()
    for i,aa in enumerate(axes):
        axes_out.append(np.matmul(R,aa))
    grp_out = list()
    grp_out.append(np.identity(3))
    for i,aa in enumerate(axes_out):
        grp_loc = generate_cyclic(aa,orders[i])
        for j,r in enumerate(grp_loc):
            if j != 0:
                grp_out.append(r)
    return(orders,axes_out,grp_out)

def tetrahedron_axes():
    axes = list()
    orders = list()
    axes.append(np.array([1,1,1.0]))
    orders.append(3)
    axes.append(np.array([1.0,1.0,-1.0]))
    orders.append(3)
    axes.append(np.array([1.0,-1.0,1.0]))
    orders.append(3)
    axes.append(np.array([-1.0,1.0,1.0]))
    orders.append(3)
    axes.append(np.array([0.0,0.0,1.0]))
    orders.append(2)
    axes.append(np.array([0.0,1.0,0.0]))
    orders.append(2)
    axes.append(np.array([1.0,0.0,0.0]))
    orders.append(2)
    for i,x in enumerate(axes):
        axes[i] = x/np.linalg.norm(x)
    return(orders,axes)

def generate_all_octahedron(axis1,order1,axis2,order2,toler=1.0e-1):
    orders,axes = octahedron_axes()
    axis_l1 = axis1/np.linalg.norm(axis1); axis_l2 = axis2/np.linalg.norm(axis2)
    if order1 == 4 and order2 == 4:
        found = 0
        for i,o in enumerate(orders):
            if o==4 and found !=2 :
                if "ax1" in locals() :
                    ax2 = axes[i]
                    found = found + 1
                else:
                    ax1 = axes[i]
                    found = found + 1
    elif order1 == 4 and order2 == 3:
        #
        # Any two four and three fold axes would do
        found = 0
        for i,o in enumerate(orders):
            if (o == 4 or o == 3)  and found != 2:
                if o == 4:
                    ax1 = axes[i]
                    found = found + 1
                else:
                    ax2 = axes[i]
                    found = found + 1
    else:
        #
        # Error
        print('Only combination of two and three or three and three folds are allowed')
        sys.exit('Wrong combinations of axes orders: must be (3,3) or (2,2)')
    ax1 = ax1/np.linalg.norm(ax1)
    ax2 = ax2/np.linalg.norm(ax2)
    #
    #  Find transformation matrix. Make sure that angles between axes is correct. Make sure that we
    #  consider ax2 and -ax2, just in case. We assume that the lengths of vectors are equal to 1.
    cangle = np.dot(axis_l1,axis_l2)
    cangle_t = np.dot(ax1,ax2)
    if np.abs(cangle+cangle_t) < toler :
        axis_l2 = -axis_l2
    ax3 = np.cross(ax1,ax2); ax3 = ax3/np.linalg.norm(ax3)
    axis_l3 = np.cross(axis_l1,axis_l2); axis_l3 = axis_l3/np.linalg.norm(axis_l3)
    AA = np.array([axis_l1,axis_l2,axis_l3])
    BB = np.array([ax1,ax2,ax3])
    from scipy.linalg import orthogonal_procrustes
    R,sca = orthogonal_procrustes(AA,BB)
    axes_out = list()
    for i,aa in enumerate(axes):
        axes_out.append(np.matmul(R,aa))
    grp_out = list()
    grp_out.append(np.identity(3))
    for i,aa in enumerate(axes_out):
        grp_loc = generate_cyclic(aa,orders[i])
        for j,r in enumerate(grp_loc):
            if j != 0:
                grp_out.append(r)
    return(orders,axes_out,grp_out)

def octahedron_axes():
    axes = list()
    orders = list()
    axes.append(np.array([0.0,0.0,1.0]))
    orders.append(4)
    axes.append(np.array([0.0,1.0,0.0]))
    orders.append(4)
    axes.append(np.array([1.0,0.0,0.0]))
    orders.append(4)
    #
    #  Three folds
    axes.append(np.array([1.0,1.0,1.0]))
    orders.append(3)
    axes.append(np.array([1.0,-1.0,1.0]))
    orders.append(3)
    axes.append(np.array([-1.0,1.0,1.0]))
    orders.append(3)
    axes.append(np.array([-1.0,-1.0,1.0]))
    orders.append(3)
    #
    # Two folds
    axes.append(np.array([1.0,0.0,1.0]))
    orders.append(2)
    axes.append(np.array([0.0,1.0,1.0]))
    orders.append(2)
    axes.append(np.array([-1.0,.0,1.0]))
    orders.append(2)
    axes.append(np.array([0.0,-1.0,0.0]))
    orders.append(2)
    axes.append(np.array([1.0,1.0,0.0]))
    orders.append(2)
    axes.append(np.array([1.0,-1.0,0.0]))
    orders.append(2)
    for i,x in enumerate(axes):
        axes[i] = x/np.linalg.norm(x)
    return(orders,axes)

def generate_all_icosahedron(axis1,order1,axis2,order2,toler=1.0e-1):
    orders,axes = icosahedron_axes()
    axis_l1 = axis1/np.linalg.norm(axis1); axis_l2 = axis2/np.linalg.norm(axis2)
    if order1 == 5 and order2 == 5:
        found = 0
        for i,o in enumerate(orders):
            if o==4 and found !=2 :
                if "ax1" in locals() :
                    ax2 = axes[i]
                    found = found + 1
                else:
                    ax1 = axes[i]
                    found = found + 1
    elif order1 == 5 and order2 == 3:
       #
       # AAngles between axes must obey certain relationship: cos(angle) = ±0.79465 or ±0.187592
       pass
    elif order1 == 5 and order2 == 2:
       #
       # Angles between axes can be: cos(anlge) = ±0.85065 or ± 0.52573
       pass
   #
   #  We can have also combinations of 3- and 2-folds as well as combinations of two 3-folds. 
   # For completeness we need to impement these also. 
#        found = 0
#        for i,o in enumerate(orders):
#                if o == 4:
#            if (o == 4 or o == 3)  and found != 2:
#                    ax1 = axes[i]
#                else:
#                    found = found + 1
#                    ax2 = axes[i]
#                    found = found + 1
    else:
        #
        # Error
        print('Only combination of five and five folds are currently allowed')
        sys.exit('Wrong combinations of axes orders: must be (5,5) ')
    ax1 = ax1/np.linalg.norm(ax1)
    ax2 = ax2/np.linalg.norm(ax2)
    #
    #  Find transformation matrix. Make sure that angles between axes is correct. Make sure that we
    #  consider ax2 and -ax2, just in case. We assume that the lengths of vectors are equal to 1.
    cangle = np.dot(axis_l1,axis_l2)
    cangle_t = np.dot(ax1,ax2)
    if np.abs(cangle+cangle_t) < toler :
        axis_l2 = -axis_l2
    ax3 = np.cross(ax1,ax2); ax3 = ax3/np.linalg.norm(ax3)
    axis_l3 = np.cross(axis_l1,axis_l2); axis_l3 = axis_l3/np.linalg.norm(axis_l3)
    AA = np.array([axis_l1,axis_l2,axis_l3])
    BB = np.array([ax1,ax2,ax3])
    from scipy.linalg import orthogonal_procrustes
    R,sca = orthogonal_procrustes(AA,BB)
    axes_out = list()
    for i,aa in enumerate(axes):
        axes_out.append(np.matmul(R,aa))
    grp_out = list()
    grp_out.append(np.identity(3))
    for i,aa in enumerate(axes_out):
        grp_loc = generate_cyclic(aa,orders[i])
        for j,r in enumerate(grp_loc):
            if j != 0:
                grp_out.append(r)
    return(orders,axes_out,grp_out)
    pass

def icosahedron_axes():
    axes = list()
    orders = list()
    verteces = list()
    #
    #  Taken from Litvin, (1990) Acta Cryst A47, pp70-73
    
    #
    # Five fold are at the vertices of the icosahedron
    u = 1.0
    r = (np.sqrt(5.0)+1.0)/2.0
    v = r-1
    verteces.append(np.array([v,0.0,u]))
    verteces.append(np.array([-v,0,u]))
    verteces.append(np.array([0.0,u,v]))
    verteces.append(np.array([u,v,0.0]))
    verteces.append(np.array([u,-v,0.0]))
    verteces.append(np.array([0.0,-u,v]))
    #
    # The rest can be derived from above multiplying by -1
    verteces.append(np.array([-u,v,0.0]))
    verteces.append(np.array([0.0,u,-v]))
    verteces.append(np.array([v,0.0,-u]))
    verteces.append(np.array([0.0,-u,-v]))
    verteces.append([-u,-v,0.0])
    verteces.append([-v,0.0,-u])
    #
    # First six are five fold axes
    for i in range(6):
        axes.append(verteces[i])
        orders.append(5)
    #
    # Three fold axes are the in the centre of faces
    ax1 = (verteces[0]+verteces[1]+verteces[2])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[0]+verteces[2]+verteces[3])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[0]+verteces[3]+verteces[4])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[0]+verteces[4]+verteces[5])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[0]+verteces[5]+verteces[1])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[2]+verteces[3]+verteces[7])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[3]+verteces[7]+verteces[8])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[3]+verteces[4]+verteces[8])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[4]+verteces[8]+verteces[9])/3.0
    axes.append(ax1)
    orders.append(3)
    ax1 = (verteces[4] + verteces[5] + verteces[9])/3.0
    axes.append(ax1)
    orders.append(3)
    #
    #  Two fold axes are on the centre of the edges.
    ax1 = (verteces[0]+verteces[1])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[0]+verteces[2])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[0]+verteces[3])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[0]+verteces[4])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[0]+verteces[5])/2.0
    axes.append(ax1)
    orders.append(2)

    ax1 = (verteces[1]+verteces[2])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[2]+verteces[3])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[3]+verteces[4])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[4]+verteces[5])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[5]+verteces[1])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[2] + verteces[7])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[3]+verteces[7])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[3]+verteces[8])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[4]+verteces[8])/2.0
    axes.append(ax1)
    orders.append(2)
    ax1 = (verteces[4]+verteces[9])/2
    axes.append(ax1)
    orders.append(2)
    for i,x in enumerate(axes):
        axes[i] = x/np.linalg.norm(x)
    return(orders,axes)

def rotate_group_elements(Rg,matrices) :
    #
    #   assume an input list and return a list of matrices
    #
    mm_out = list()
    Rg_t = np.transpose(Rg)
    for i,mm in enumerate(matrices) :
        m1 = np.dot(np.dot(Rg_t,mm),Rg)
        mm_out.append(m1)
    return(mm_out)
    
def get_matrices_using_relion(sym):
    import subprocess
    ps = subprocess.check_output(["relion_refine", "--sym", sym.strip(), "--print_symmetry_ops"])

    ret = []
    read_flag = -1
    for l in ps.splitlines():
        if b"R(" in l:
            ret.append(np.zeros((3,3)))
            read_flag = 0
        elif 0 <= read_flag < 3:
            ret[-1][read_flag,:] = [float(x) for x in l.split()]
            read_flag += 1
        elif read_flag >= 3:
            read_flag = -1

    return ret
# get_matrices_using_relion()

if __name__ == "__main__":
    import sys
    symbol = sys.argv[1]
    order, axes, grp = operators_from_symbol(symbol)
    print(order)
    print(axes)
    print(grp)
    
    '''rgs = get_matrices_using_relion(symbol)
    all_ok = True
    max_diff = 0
    for i, m in enumerate(grp):
        print("Op", i)
        print(m)
        #ok = is_in_the_list_rotation(m, rgs, toler=1e-4)
        diff = closest_rotation(m, rgs)
        ok = diff < 1e-4
        print("match? {} {:.1e}".format(ok, diff))
        if not ok: all_ok = False
        if diff > max_diff: max_diff = diff
    print("Final=", all_ok, max_diff)'''
