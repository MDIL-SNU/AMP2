###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
# This is a package of modules to calculate several properties such as distance between two points and angle.

# Return reciprocal vectors from POSCAR
def reciprocal_lattice(axis) :
	bm = (axis[1][1]*axis[2][2]-axis[1][2]*axis[2][1])*axis[0][0]+(axis[1][2]*axis[2][0]-axis[1][0]*axis[2][2])*axis[0][1]+(axis[1][0]*axis[2][1]-axis[1][1]*axis[2][0])*axis[0][2]
	b1 = [(axis[1][1]*axis[2][2]-axis[1][2]*axis[2][1])/bm, (axis[1][2]*axis[2][0]-axis[1][0]*axis[2][2])/bm, (axis[1][0]*axis[2][1]-axis[1][1]*axis[2][0])/bm]
	bm = -(axis[0][1]*axis[2][2]-axis[0][2]*axis[2][1])*axis[1][0]-(axis[0][2]*axis[2][0]-axis[0][0]*axis[2][2])*axis[1][1]-(axis[0][0]*axis[2][1]-axis[0][1]*axis[2][0])*axis[1][2]
	b2 = [-(axis[0][1]*axis[2][2]-axis[0][2]*axis[2][1])/bm, -(axis[0][2]*axis[2][0]-axis[0][0]*axis[2][2])/bm, -(axis[0][0]*axis[2][1]-axis[0][1]*axis[2][0])/bm]
	bm = (axis[0][1]*axis[1][2]-axis[0][2]*axis[1][1])*axis[2][0]+(axis[0][2]*axis[1][0]-axis[0][0]*axis[1][2])*axis[2][1]+(axis[0][0]*axis[1][1]-axis[0][1]*axis[1][0])*axis[2][2]
	b3 = [(axis[0][1]*axis[1][2]-axis[0][2]*axis[1][1])/bm, (axis[0][2]*axis[1][0]-axis[0][0]*axis[1][2])/bm, (axis[0][0]*axis[1][1]-axis[0][1]*axis[1][0])/bm]

	return [b1,b2,b3]

# Calculate distance between two vector
def dist_vec(v1,v2,lat) :
	vd = [float(x)-float(y) for x,y in zip(v1,v2)]
	dist = ((lat[0][0]*vd[0]+lat[1][0]*vd[1]+lat[2][0]*vd[2])**2.0+(lat[0][1]*vd[0]+lat[1][1]*vd[1]+lat[2][1]*vd[2])**2.0+(lat[0][2]*vd[0]+lat[1][2]*vd[1]+lat[2][2]*vd[2])**2.0)**0.5
	return dist

# This function is for calculating distance between two points
def dist_point(p1,p2) :
	dist = 0
	for i in range(len(p1)):
		dist = dist+(float(p1[i])-float(p2[i]))**2.0
	dist = dist**0.5
	return dist

# This function is for converting from direct coordinate to cartesian coordinate
def dir_to_cart(vec,axis):
	cart_vec = []
	for i in range(3):
		cart_vec.append(vec[0]*axis[0][i]+vec[1]*axis[1][i]+vec[2]*axis[2][i])
	return cart_vec

# This function is for converting from cartesian coordinate to direct coordinate
def cart_to_dir(vec,axis):
	import numpy as np
	vec_array = np.array(vec)
	axis_array = np.array(axis)
	dir_vec_array = np.dot(vec_array,np.linalg.inv(axis_array))
	dir_vec = dir_vec_array.tolist()
	return dir_vec

# This function is for calculating angle between two vectors
def calc_angle(vec1,vec2):
	from math import acos
	dot = sum([float(vec1[x])*float(vec2[x]) for x in range(len(vec1))])
	len1 = dist_point(vec1,[0,0,0])
	len2 = dist_point(vec2,[0,0,0])
	costheta = dot/(len1*len2)
	theta = acos(costheta)
	return theta

# This function is for calculating distance applying periodic boundary condition
def short_dist(v1,v2,lat) :
	dup_list=[-2,-1,0,1,2]
	short_length = 999
	for x in dup_list:
		for y in dup_list:
			for z in dup_list:
				dup = [x,y,z]
				length = dist_vec([float(v1[k]) + float(dup[k]) for k  in range(3)], v2,lat)
				if short_length > length and length > 0.1:
					short_length = length
	return short_length

# This function is for calculating volume of structure
def calc_volume(axis):
	import numpy as np
	axis = np.array(axis)
	volume = np.dot(axis[0],np.cross(axis[1],axis[2]))
	return volume
