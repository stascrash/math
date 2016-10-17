import math


class Vec(object):
	def __init__(self, x=0, y=0, z=0):
		self.x = x
		self.y = y
		self.z = z

	def __call__(self, *args, **kwargs):
		return self.x, self.y, self.z


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



class Circle(object):
	def __init__(self, r, x, y):
		self.center = Point(x, y)
		self.radius = r

	def intersect(self, c2):
		# c1 - Circle One, c2 - Circule Two

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

			chord_center = self.center + (c2.center - self.center)*a / distance

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

		self.points = [self.A, self.B, self.C] # adding to vertex-list for convenience
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

		# Distance from each vertex to centroid
		# self.a_centroid, self.b_centroid, self.c_centroid = self.get_centroid_lengths()

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
		Ox = (self.A.x + self.B.x + self.C.x) / 3.
		Oy = (self.A.y + self.B.y + self.C.y) / 3.
		return Point(x=Ox, y=Oy)

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
		circle_a = Circle(side_b, pt_a.x, pt_a.y)

		# Plot point B on X-axis with distance x=side_a
		pt_b = Point(x=side_a)
		# Build circle with center at point B and radius side_c
		circle_b = Circle(side_c, pt_b.x, pt_b.y)

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


def translate_coords(object_, new_coords_point):
	"""
	Change coordinate for object to new. adjust all points to
	represent their position in new coordinate system.
	"""
	for point in object_.points:
		point.x = point.x - new_coords_point.x
		point.y = point.y - new_coords_point.y
		point.z = point.z - new_coords_point.z

def rotate_obj_on_xyz(object_, x_axis_deg, y_axis_deg, z_axis_deg):
	for point in object_.points:
		# Find current direction relative to center
		distance = point.get_distance(object_.centroid)

		# A, B, G - alpha, Beta, Gamma angle abbreviations
		# Get direction-angle for each axis
		cosA = point.x / distance
		cosB = point.y / distance
		cosG = point.z / distance
		# Set global rotation in new position (radians) and get cos
		newCosA = math.cos(math.acos(cosA) + math.radians(x_axis_deg))
		newCosB = math.cos(math.acos(cosB) + math.radians(y_axis_deg))
		newCosG = math.cos(math.acos(cosG) + math.radians(z_axis_deg))

		# Set new X,Y,Z for point after rotation
		point.x = distance * newCosA
		point.y = distance * newCosB
		point.z = distance * newCosG

def main():
	triangle = Triangle(6, 6, 6)
	new_sys = triangle.centroid
	new_sys.z = 1.
	translate_coords(triangle, triangle.centroid)
	triangle.update_properties()
	rotate_obj_on_xyz(triangle, 60, 30, 30)
	triangle.update_properties()
	translate_coords(triangle, -new_sys)
	triangle.update_properties()
	print(triangle)



if __name__ == '__main__':
	main()





