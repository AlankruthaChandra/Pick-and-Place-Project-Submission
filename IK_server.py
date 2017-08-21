#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np
from numpy import array
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos
from sympy.matrices import Matrix


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
	joint_trajectory_list = []
    
		
        d1,d2,d3,d4,d5,d6,d7 = symbols('d1:8')
    	a0,a1,a2,a3,a4,a5,a6 = symbols('a0:7')
    	alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')


    # Joint angle symbols
    	q1,q2,q3,q4,q5,q6,q7 = symbols('q1:8')


    # Modified DH params
    	s={     alpha0:		0,a0:           0, d1:          0.75,q1:  q1,
   		alpha1:     -pi/2,a1:	     0.35, d2:		   0,q2: -pi/2+q2, 
   		alpha2:		0,a2:	     1.25, d3:		   0,       
   		alpha3:     -pi/2,a3:      -0.054, d4:          1.50,       
   		alpha4:      pi/2,a4:          0,  d5:             0,         
   		alpha5:     -pi/2,a5:          0,  d6:             0,         
   		alpha6:        0, a6:          0,  d7:         0.303,q7:          0}

    # Define Modified DH Transformation matrix
        def TF_Matrix(alpha, a, d, q):
	    TF = Matrix([[     	     cos(q),     	 -sin(q),            0,              a],
                     [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
              	     [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
               	     [                   0,                   0,            0,               1]])
	    return TF

	    
    #Create individual transformation matrices
        T0_1 = TF_Matrix(alpha0, a0, d1, q1).subs(s)
    	T1_2 = TF_Matrix(alpha1, a1, d2, q2).subs(s)
    	T2_3 = TF_Matrix(alpha2, a2, d3, q3).subs(s)
    	T3_4 = TF_Matrix(alpha3, a3, d4, q4).subs(s)
    	T4_5 = TF_Matrix(alpha4, a4, d5, q5).subs(s)
    	T5_6 = TF_Matrix(alpha5, a5, d6, q6).subs(s)
    	T6_7 = TF_Matrix(alpha6, a6, d7, q7).subs(s)

    	    
    	#T0_3 = simplify(T0_1 * T1_2 * T2_3)
    	    #T0_6 = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 
	    #T0_7 = T0_6 * T6_7


    # Extract end-effector position and orientation from request
    # px,py,pz = end-effector position
    # roll, pitch, yaw = end-effector orientation
    	for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()
	    p_x = req.poses[x].position.x
    	    p_y = req.poses[x].position.y
            p_z = req.poses[x].position.z

    	    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            # Calculate joint angles using Geometric IK method
    	    r, p, y = symbols('r p y')
    
            R_x  = Matrix([[1,      0,       0],
		   [0, cos(r), -sin(r)],
		   [0 ,sin(r),  cos(r)]])

            R_y  = Matrix([[cos(p), 0, sin(p)],
		   [0,      1,       0],
		   [-sin(p), 0,  cos(p)]])

            R_z  = Matrix([[cos(y), -sin(y), 0],
		   [sin(y),  cos(y), 0],
		   [0 ,     0,       1]])
            Rrpy = R_z * R_y * R_x

    	    R_correction = R_z.subs(y, pi) * R_y.subs(p, -pi/2)


            Rrpy = Rrpy. subs({'r': roll, 'p' :pitch, 'y': yaw})

            Rrpy = Rrpy* R_correction

            R = Matrix([[p_x],
	        [p_y],
	        [p_z]])
   	    #wrist center calculation
            WC = (R - (0.303) * Matrix([[Rrpy[0,2]],
			       [Rrpy[1,2]],
			       [Rrpy[2,2]]])).evalf()
    

    	    a = 1.501 #s[d4]
            b = sqrt(pow((sqrt(WC[0]*WC[0] + WC[1] * WC[1])-0.35),2)+pow((WC[2]-0.75),2))
            c = 1.25 #s[a2]

	    
    
            theta1_2 = acos((b*b + c*c - a*a)/(2*b*c))
            theta2_3 = acos((a*a + c*c - b*b)/(2*a*c))
  
            theta1 = atan2(WC[1],WC[0]).evalf()
            theta2 = (pi/2 -  theta1_2 - atan2(WC[2] - 0.75, sqrt(WC[0]*WC[0] + WC[1] * WC[1]) - 0.35)).evalf()
            theta3 = (pi/2 - (theta2_3 + 0.036)).evalf()
	   
            R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
            R0_3 = R0_3.evalf(subs={q1:theta1, q2:theta2, q3:theta3})
	    
            R3_6 = R0_3.inv("LU") * Rrpy

	    
	    theta5 = atan2(sqrt(pow(R3_6[0,2],2)+pow(R3_6[2,2],2)), R3_6[1,2])
	    #theta5 = atan2(sqrt(pow(R3_6[1,0],2)+pow(R3_6[1,1],2)), R3_6[1,2])
	    if sin(theta5) < 0:
		theta4 = atan2(R3_6[2,2], -R3_6[0,2])
	        theta6 = atan2(-R3_6[1,1], R3_6[1,0])

	    else:
		theta4 = atan2(-R3_6[2,2], R3_6[0,2])
		theta6 = atan2(R3_6[1,1], -R3_6[1,0])
	   
	    
	    #theta5 = atan2(sqrt(pow(R3_6[1,0],2)+pow(R3_6[1,1],2)), R3_6[1,2])
	    #theta4 = atan2(sqrt(R))
	    #theta6 = atan2(-R3_6[1,1], R3_6[1,0]) 


	
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
    	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
    	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
