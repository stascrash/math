import math


def translate_coords(object_, new_coords_point):
	"""
	Change coordinate for object to new. adjust all points to
	represent their position in new coordinate system.
	"""
	for point in object_.points:
		point.x = point.x - new_coords_point.x
		point.y = point.y - new_coords_point.y
		point.z = point.z - new_coords_point.z


def rotate_obj_on_xyz(object_, rot_vec):
	for point in object_.points:
		# Find current direction relative to center
		distance = point.get_distance(object_.centroid)

		# A, B, G - alpha, Beta, Gamma angle abbreviations
		# Get direction-angle for each axis
		cos_a = point.x / distance
		cos_b = point.y / distance
		cos_g = point.z / distance
		# Set global rotation in new position (radians) and get cos
		new_cos_a = math.cos(math.acos(cos_a) + math.radians(rot_vec.x))
		new_cos_b = math.cos(math.acos(cos_b) + math.radians(rot_vec.y))
		new_cos_g = math.cos(math.acos(cos_g) + math.radians(rot_vec.z))

		# Set new X,Y,Z for point after rotation
		point.x = distance * new_cos_a
		point.y = distance * new_cos_b
		point.z = distance * new_cos_g


class Point(object):
	def __init__(self, x=0., y=0., z=0.):
		self.x = x
		self.y = y
		self.z = z

	def get_distance(self, B):
		return math.sqrt((B.x - self.x)**2 + (B.y - self.y)**2 + (B.z - self.z)**2)

	def __add__(self, other):
		x = self.x + other.x
		y = self.y + other.y
		z = self.z + other.z
		return Point(x, y, z)

	def __sub__(self, other):
		x = self.x - other.x
		y = self.y - other.y
		z = self.z - other.z
		return Point(x, y, z)

	def __mul__(self, other):
		x = self.x * other
		y = self.y * other
		z = self.z * other
		return Point(x, y, z)

	def __floordiv__(self, other):
		x = self.x // other
		y = self.y // other
		z = self.z // other
		return Point(x, y, z)

	def __neg__(self):
		return Point(-self.x, -self.y, -self.z)

	def __truediv__(self, other):
		x = self.x / other
		y = self.y / other
		z = self.z / other
		return Point(x, y, z)


class Vec(Point):
	def __init__(self, x=0, y=0, z=0):
		super(Vec, self).__init__(x, y, z)

	def __call__(self, *args, **kwargs):
		return self.x, self.y, self.z


class Circle(object):
	def __init__(self, r, x, y, z):
		self.center = Point(x, y, z)
		self.radius = r

	def intersect(self, c2):
		# self - Circle One, c2 - Circle Two

		# Get distance between circle-centers
		distance = math.sqrt((c2.center.x - self.center.x) ** 2 + (c2.center.y - self.center.y) ** 2)

		# Check if circles can be intersected at all
		if distance > (self.radius + c2.radius):
			raise Exception('There are no solutions, the circles are separate')
		elif distance < abs((self.radius - c2.radius)):
			raise Exception('There are no solutions because one circle is contained within the other')
		elif distance == 0:
			raise Exception('The circles are coincident and there are an infinite number of solutions')
		else:
			a = (self.radius**2 - c2.radius**2 + distance**2) / (2 * distance)
			b = distance - a
			h = math.sqrt(self.radius**2 - a**2)

			chord_center = self.center + (c2.center - self.center) * a / distance

			intersection_p1_x = chord_center.x + h * (self.center.y - c2.center.y) / distance
			intersection_p1_y = chord_center.y - h * (self.center.x - c2.center.x) / distance
			p1 = Point(x=intersection_p1_x, y=intersection_p1_y)

			intersection_p2_x = chord_center.x - h * (self.center.y - c2.center.y) / distance
			intersection_p2_y = chord_center.y + h * (self.center.x - c2.center.x) / distance
			p2 = Point(x=intersection_p2_x, y=intersection_p2_y)

			return p1, p2


class Triangle(object):
	"""
	Triangle object by length of sides.
	Builds a triangle with A, B, C vertexes.
	Starting vertex is A and is located in XOY plane with A(0,0,0),
	where side AB is a base and belongs to X-Axis.
	"""

	def __init__(self, *size):
		# Vertexes:
		self.A, self.B, self.C = self.get_vertexes(*size)

		self.points = [self.A, self.B, self.C]  # adding to vertex-list for convenience
		self.update_properties()

	def update_properties(self):
		"""Re-calculate all properties"""
		# Sides Distance
		self.AB, self.AC, self.BC = self.get_sides_length()
		# Angles in degrees
		self.alpha, self.beta, self.gamma = self.get_angles()

		self.centroid = self.get_center()
		self.perimeter = self.get_perimeter()
		self.half_perimeter = self.perimeter / 2.
		self.area = self.get_area()
		self.height = self.get_height()
		self.height_pt = self.get_height_point()

	def get_centroid_lengths(self):
		"""return distance from each vertex to centroid"""
		a_center = self.centroid.get_distance(self.A)
		b_center = self.centroid.get_distance(self.B)
		c_center = self.centroid.get_distance(self.C)
		return a_center, b_center, c_center

	def get_height_point(self):
		AD = self.height # Length of side AD
		BD = math.sqrt((self.AB**2 - AD**2))
		CD = abs(self.BC - BD)
		ratio = BD / CD
		x = (self.B.x + ratio*self.C.x) / (1 + ratio)
		y = (self.B.y + ratio*self.C.y) / (1 + ratio)
		return Point(x, y)

	def get_height(self):
		return (2 * self.area) / self.BC

	def get_area(self):
		return math.sqrt((self.half_perimeter*((self.half_perimeter - self.AB)*(self.half_perimeter - self.BC)*(self.half_perimeter - self.AC))))

	def get_perimeter(self):
		return self.AB + self.BC + self.AC

	def get_center(self):
		# Center of side CB Where D is mid-point, Xd - coord for X and Yd coord for Y
		# Center of Triangle: P, Xp - for coord X and Yp for coord Y
		# Distance from A to D, therefore lambda is 2
		Xd = (self.B.x + self.C.x) / 2.
		Yd = (self.B.y + self.C.y) / 2.
		Zd = (self.B.z + self.C.z) / 2.
		Xp = (self.A.x + 2 * Xd) / 3.
		Yp = (self.A.y + 2 * Yd) / 3.
		Zp = (self.A.z + 2 * Zd) / 3.
		return Point(Xp, Yp, Zp)

	def get_vertexes(self, *size):
		"""
		Build a triangle project on plane: XOY
		Where: x - Axis X | y - Axis Y | O - Start of coordinate plane.
		"""
		side_a, side_b, side_c = size
		# Check if supplied sizes can be used to build triangle:
		if (side_a + side_b) < side_c or \
			(side_b + side_c) < side_a or \
			(side_a + side_c) < side_b:
			raise ValueError('Given sizes cannot form triangle')

		# Plot  point A at origin 0,0
		pt_a = Point()
		# Build circle - A with center in pt_a and radius of side_b
		circle_a = Circle(side_b, pt_a.x, pt_a.y, pt_a.z)

		# Plot point B on X-axis with distance x=side_a
		pt_b = Point(x=side_a)
		# Build circle with center at point B and radius side_c
		circle_b = Circle(side_c, pt_b.x, pt_b.y, pt_b.z)

		# Intersect Circle_A with Circle_B to get (x,y) location for pt_c
		pt_c = circle_a.intersect(circle_b)[0]

		return pt_a, pt_b, pt_c

	def get_sides_length(self):
		ac = self.A.get_distance(self.C)
		bc = self.B.get_distance(self.C)
		ab = self.A.get_distance(self.B)
		return ab, ac, bc

	def get_angles(self):
		alpha = (self.AC**2 + self.AB**2 - self.BC**2) / (2 * self.AC * self.AB)
		beta = (self.AB**2 + self.BC**2 - self.AC**2) / (2 * self.AB * self.BC)
		gamma = (self.BC**2 + self.AC**2 - self.AB**2) / (2 * self.BC * self.AC)
		return math.degrees(math.acos(alpha)), math.degrees(math.acos(beta)), math.degrees(math.acos(gamma))


class Prism(object):
	def __init__(self, height, *size):
		self._top = Triangle(*size)
		translate_coords(self._top, Point(0, 0, height))
		self._bottom = Triangle(*size)
		self.points = list()
		self.points.extend(self._top.points)
		self.points.extend(self._bottom.points)
		self.centroid = self._top.centroid
		self.centroid.z = height / 2.

	def update(self):
		self._top.update_properties()
		self._bottom.update_properties()


def rotate_prism(rotation_vec):
	prism = Prism(2, 6, 6, 6)
	coord_origin = prism.centroid

	translate_coords(prism, coord_origin)
	prism.update()

	rotate_obj_on_xyz(prism, rotation_vec)
	prism.update()

	translate_coords(prism, -coord_origin)
	prism.update()
	return prism


def rotate_triangles(rotation_vec):
	# Create triangle pair
	tri_a = Triangle(6, 6, 6)
	# tri_b = Triangle(6, 6, 6)

	# Translate 2 units up
	# translate_coords(tri_b, Point(0, 0, 2))
	# tri_b.update_properties()

	# Get point representing new origin of coordinate system
	coord_origin = tri_a.centroid
	# coord_origin.z = 1.

	# translate all triangle-point to work with new origin
	translate_coords(tri_a, coord_origin)
	# translate_coords(tri_b, coord_origin)
	tri_a.update_properties()
	# tri_b.update_properties()

	# Rotate Triangle's points
	rotate_obj_on_xyz(tri_a, rotation_vec)
	# rotate_obj_on_xyz(tri_b, rotation_vec)
	tri_a.update_properties()
	# tri_b.update_properties()

	# Translate triangle back to original coordinate system
	translate_coords(tri_a, -coord_origin)
	# translate_coords(tri_b, -coord_origin)
	tri_a.update_properties()
	# tri_b.update_properties()
	return tri_a  #, tri_b


def main():
	rotation_vec = Vec(60, 30, 30)
	tri_a= rotate_triangles(rotation_vec)
	print('Finish')




if __name__ == '__main__':
	main()





