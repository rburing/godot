/*************************************************************************/
/*  geometry_3d.cpp                                                      */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2022 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2022 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "geometry_3d.h"

#include "polynomial.h" // for circle-circle
#include "polynomial_root_finder.h" // for line-circle and circle-circle

#include "thirdparty/misc/clipper.hpp"
#include "thirdparty/misc/polypartition.h"

void Geometry3D::get_closest_points_between_segments(const Vector3 &p_p0, const Vector3 &p_p1, const Vector3 &p_q0, const Vector3 &p_q1, Vector3 &r_ps, Vector3 &r_qt) {
	// Based on David Eberly's Computation of Distance Between Line Segments algorithm.

	Vector3 p = p_p1 - p_p0;
	Vector3 q = p_q1 - p_q0;
	Vector3 r = p_p0 - p_q0;

	real_t a = p.dot(p);
	real_t b = p.dot(q);
	real_t c = q.dot(q);
	real_t d = p.dot(r);
	real_t e = q.dot(r);

	real_t s = 0.0f;
	real_t t = 0.0f;

	real_t det = a * c - b * b;
	if (det > CMP_EPSILON) {
		// Non-parallel segments
		real_t bte = b * e;
		real_t ctd = c * d;

		if (bte <= ctd) {
			// s <= 0.0f
			if (e <= 0.0f) {
				// t <= 0.0f
				s = (-d >= a ? 1 : (-d > 0.0f ? -d / a : 0.0f));
				t = 0.0f;
			} else if (e < c) {
				// 0.0f < t < 1
				s = 0.0f;
				t = e / c;
			} else {
				// t >= 1
				s = (b - d >= a ? 1 : (b - d > 0.0f ? (b - d) / a : 0.0f));
				t = 1;
			}
		} else {
			// s > 0.0f
			s = bte - ctd;
			if (s >= det) {
				// s >= 1
				if (b + e <= 0.0f) {
					// t <= 0.0f
					s = (-d <= 0.0f ? 0.0f : (-d < a ? -d / a : 1));
					t = 0.0f;
				} else if (b + e < c) {
					// 0.0f < t < 1
					s = 1;
					t = (b + e) / c;
				} else {
					// t >= 1
					s = (b - d <= 0.0f ? 0.0f : (b - d < a ? (b - d) / a : 1));
					t = 1;
				}
			} else {
				// 0.0f < s < 1
				real_t ate = a * e;
				real_t btd = b * d;

				if (ate <= btd) {
					// t <= 0.0f
					s = (-d <= 0.0f ? 0.0f : (-d >= a ? 1 : -d / a));
					t = 0.0f;
				} else {
					// t > 0.0f
					t = ate - btd;
					if (t >= det) {
						// t >= 1
						s = (b - d <= 0.0f ? 0.0f : (b - d >= a ? 1 : (b - d) / a));
						t = 1;
					} else {
						// 0.0f < t < 1
						s /= det;
						t /= det;
					}
				}
			}
		}
	} else {
		// Parallel segments
		if (e <= 0.0f) {
			s = (-d <= 0.0f ? 0.0f : (-d >= a ? 1 : -d / a));
			t = 0.0f;
		} else if (e >= c) {
			s = (b - d <= 0.0f ? 0.0f : (b - d >= a ? 1 : (b - d) / a));
			t = 1;
		} else {
			s = 0.0f;
			t = e / c;
		}
	}

	r_ps = (1 - s) * p_p0 + s * p_p1;
	r_qt = (1 - t) * p_q0 + t * p_q1;
}

real_t Geometry3D::get_closest_distance_between_segments(const Vector3 &p_p0, const Vector3 &p_p1, const Vector3 &p_q0, const Vector3 &p_q1) {
	Vector3 ps;
	Vector3 qt;
	get_closest_points_between_segments(p_p0, p_p1, p_q0, p_q1, ps, qt);
	Vector3 st = qt - ps;
	return st.length();
}

void Geometry3D::MeshData::optimize_vertices() {
	HashMap<int, int> vtx_remap;

	for (uint32_t i = 0; i < faces.size(); i++) {
		for (uint32_t j = 0; j < faces[i].indices.size(); j++) {
			int idx = faces[i].indices[j];
			if (!vtx_remap.has(idx)) {
				int ni = vtx_remap.size();
				vtx_remap[idx] = ni;
			}

			faces[i].indices[j] = vtx_remap[idx];
		}
	}

	for (uint32_t i = 0; i < edges.size(); i++) {
		int a = edges[i].vertex_a;
		int b = edges[i].vertex_b;

		if (!vtx_remap.has(a)) {
			int ni = vtx_remap.size();
			vtx_remap[a] = ni;
		}
		if (!vtx_remap.has(b)) {
			int ni = vtx_remap.size();
			vtx_remap[b] = ni;
		}

		edges[i].vertex_a = vtx_remap[a];
		edges[i].vertex_b = vtx_remap[b];
	}

	LocalVector<Vector3> new_vertices;
	new_vertices.resize(vtx_remap.size());

	for (uint32_t i = 0; i < vertices.size(); i++) {
		if (vtx_remap.has(i)) {
			new_vertices[vtx_remap[i]] = vertices[i];
		}
	}
	vertices = new_vertices;
}

struct _FaceClassify {
	struct _Link {
		int face = -1;
		int edge = -1;
		void clear() {
			face = -1;
			edge = -1;
		}
		_Link() {}
	};
	bool valid = false;
	int group = -1;
	_Link links[3];
	Face3 face;
	_FaceClassify() {}
};

static bool _connect_faces(_FaceClassify *p_faces, int len, int p_group) {
	// Connect faces, error will occur if an edge is shared between more than 2 faces.
	// Clear connections.

	bool error = false;

	for (int i = 0; i < len; i++) {
		for (int j = 0; j < 3; j++) {
			p_faces[i].links[j].clear();
		}
	}

	for (int i = 0; i < len; i++) {
		if (p_faces[i].group != p_group) {
			continue;
		}
		for (int j = i + 1; j < len; j++) {
			if (p_faces[j].group != p_group) {
				continue;
			}

			for (int k = 0; k < 3; k++) {
				Vector3 vi1 = p_faces[i].face.vertex[k];
				Vector3 vi2 = p_faces[i].face.vertex[(k + 1) % 3];

				for (int l = 0; l < 3; l++) {
					Vector3 vj2 = p_faces[j].face.vertex[l];
					Vector3 vj1 = p_faces[j].face.vertex[(l + 1) % 3];

					if (vi1.distance_to(vj1) < 0.00001f &&
							vi2.distance_to(vj2) < 0.00001f) {
						if (p_faces[i].links[k].face != -1) {
							ERR_PRINT("already linked\n");
							error = true;
							break;
						}
						if (p_faces[j].links[l].face != -1) {
							ERR_PRINT("already linked\n");
							error = true;
							break;
						}

						p_faces[i].links[k].face = j;
						p_faces[i].links[k].edge = l;
						p_faces[j].links[l].face = i;
						p_faces[j].links[l].edge = k;
					}
				}
				if (error) {
					break;
				}
			}
			if (error) {
				break;
			}
		}
		if (error) {
			break;
		}
	}

	for (int i = 0; i < len; i++) {
		p_faces[i].valid = true;
		for (int j = 0; j < 3; j++) {
			if (p_faces[i].links[j].face == -1) {
				p_faces[i].valid = false;
			}
		}
	}
	return error;
}

static bool _group_face(_FaceClassify *p_faces, int len, int p_index, int p_group) {
	if (p_faces[p_index].group >= 0) {
		return false;
	}

	p_faces[p_index].group = p_group;

	for (int i = 0; i < 3; i++) {
		ERR_FAIL_INDEX_V(p_faces[p_index].links[i].face, len, true);
		_group_face(p_faces, len, p_faces[p_index].links[i].face, p_group);
	}

	return true;
}

Vector<Vector<Face3>> Geometry3D::separate_objects(Vector<Face3> p_array) {
	Vector<Vector<Face3>> objects;

	int len = p_array.size();

	const Face3 *arrayptr = p_array.ptr();

	Vector<_FaceClassify> fc;

	fc.resize(len);

	_FaceClassify *_fcptr = fc.ptrw();

	for (int i = 0; i < len; i++) {
		_fcptr[i].face = arrayptr[i];
	}

	bool error = _connect_faces(_fcptr, len, -1);

	ERR_FAIL_COND_V_MSG(error, Vector<Vector<Face3>>(), "Invalid geometry.");

	// Group connected faces in separate objects.

	int group = 0;
	for (int i = 0; i < len; i++) {
		if (!_fcptr[i].valid) {
			continue;
		}
		if (_group_face(_fcptr, len, i, group)) {
			group++;
		}
	}

	// Group connected faces in separate objects.

	for (int i = 0; i < len; i++) {
		_fcptr[i].face = arrayptr[i];
	}

	if (group >= 0) {
		objects.resize(group);
		Vector<Face3> *group_faces = objects.ptrw();

		for (int i = 0; i < len; i++) {
			if (!_fcptr[i].valid) {
				continue;
			}
			if (_fcptr[i].group >= 0 && _fcptr[i].group < group) {
				group_faces[_fcptr[i].group].push_back(_fcptr[i].face);
			}
		}
	}

	return objects;
}

/*** GEOMETRY WRAPPER ***/

enum _CellFlags {
	_CELL_SOLID = 1,
	_CELL_EXTERIOR = 2,
	_CELL_STEP_MASK = 0x1C,
	_CELL_STEP_NONE = 0 << 2,
	_CELL_STEP_Y_POS = 1 << 2,
	_CELL_STEP_Y_NEG = 2 << 2,
	_CELL_STEP_X_POS = 3 << 2,
	_CELL_STEP_X_NEG = 4 << 2,
	_CELL_STEP_Z_POS = 5 << 2,
	_CELL_STEP_Z_NEG = 6 << 2,
	_CELL_STEP_DONE = 7 << 2,
	_CELL_PREV_MASK = 0xE0,
	_CELL_PREV_NONE = 0 << 5,
	_CELL_PREV_Y_POS = 1 << 5,
	_CELL_PREV_Y_NEG = 2 << 5,
	_CELL_PREV_X_POS = 3 << 5,
	_CELL_PREV_X_NEG = 4 << 5,
	_CELL_PREV_Z_POS = 5 << 5,
	_CELL_PREV_Z_NEG = 6 << 5,
	_CELL_PREV_FIRST = 7 << 5,
};

static inline void _plot_face(uint8_t ***p_cell_status, int x, int y, int z, int len_x, int len_y, int len_z, const Vector3 &voxelsize, const Face3 &p_face) {
	AABB aabb(Vector3(x, y, z), Vector3(len_x, len_y, len_z));
	aabb.position = aabb.position * voxelsize;
	aabb.size = aabb.size * voxelsize;

	if (!p_face.intersects_aabb(aabb)) {
		return;
	}

	if (len_x == 1 && len_y == 1 && len_z == 1) {
		p_cell_status[x][y][z] = _CELL_SOLID;
		return;
	}

	int div_x = len_x > 1 ? 2 : 1;
	int div_y = len_y > 1 ? 2 : 1;
	int div_z = len_z > 1 ? 2 : 1;

#define SPLIT_DIV(m_i, m_div, m_v, m_len_v, m_new_v, m_new_len_v) \
	if (m_div == 1) {                                             \
		m_new_v = m_v;                                            \
		m_new_len_v = 1;                                          \
	} else if (m_i == 0) {                                        \
		m_new_v = m_v;                                            \
		m_new_len_v = m_len_v / 2;                                \
	} else {                                                      \
		m_new_v = m_v + m_len_v / 2;                              \
		m_new_len_v = m_len_v - m_len_v / 2;                      \
	}

	int new_x;
	int new_len_x;
	int new_y;
	int new_len_y;
	int new_z;
	int new_len_z;

	for (int i = 0; i < div_x; i++) {
		SPLIT_DIV(i, div_x, x, len_x, new_x, new_len_x);

		for (int j = 0; j < div_y; j++) {
			SPLIT_DIV(j, div_y, y, len_y, new_y, new_len_y);

			for (int k = 0; k < div_z; k++) {
				SPLIT_DIV(k, div_z, z, len_z, new_z, new_len_z);

				_plot_face(p_cell_status, new_x, new_y, new_z, new_len_x, new_len_y, new_len_z, voxelsize, p_face);
			}
		}
	}

#undef SPLIT_DIV
}

static inline void _mark_outside(uint8_t ***p_cell_status, int x, int y, int z, int len_x, int len_y, int len_z) {
	if (p_cell_status[x][y][z] & 3) {
		return; // Nothing to do, already used and/or visited.
	}

	p_cell_status[x][y][z] = _CELL_PREV_FIRST;

	while (true) {
		uint8_t &c = p_cell_status[x][y][z];

		if ((c & _CELL_STEP_MASK) == _CELL_STEP_NONE) {
			// Haven't been in here, mark as outside.
			p_cell_status[x][y][z] |= _CELL_EXTERIOR;
		}

		if ((c & _CELL_STEP_MASK) != _CELL_STEP_DONE) {
			// If not done, increase step.
			c += 1 << 2;
		}

		if ((c & _CELL_STEP_MASK) == _CELL_STEP_DONE) {
			// Go back.
			switch (c & _CELL_PREV_MASK) {
				case _CELL_PREV_FIRST: {
					return;
				} break;
				case _CELL_PREV_Y_POS: {
					y++;
					ERR_FAIL_COND(y >= len_y);
				} break;
				case _CELL_PREV_Y_NEG: {
					y--;
					ERR_FAIL_COND(y < 0);
				} break;
				case _CELL_PREV_X_POS: {
					x++;
					ERR_FAIL_COND(x >= len_x);
				} break;
				case _CELL_PREV_X_NEG: {
					x--;
					ERR_FAIL_COND(x < 0);
				} break;
				case _CELL_PREV_Z_POS: {
					z++;
					ERR_FAIL_COND(z >= len_z);
				} break;
				case _CELL_PREV_Z_NEG: {
					z--;
					ERR_FAIL_COND(z < 0);
				} break;
				default: {
					ERR_FAIL();
				}
			}
			continue;
		}

		int next_x = x, next_y = y, next_z = z;
		uint8_t prev = 0;

		switch (c & _CELL_STEP_MASK) {
			case _CELL_STEP_Y_POS: {
				next_y++;
				prev = _CELL_PREV_Y_NEG;
			} break;
			case _CELL_STEP_Y_NEG: {
				next_y--;
				prev = _CELL_PREV_Y_POS;
			} break;
			case _CELL_STEP_X_POS: {
				next_x++;
				prev = _CELL_PREV_X_NEG;
			} break;
			case _CELL_STEP_X_NEG: {
				next_x--;
				prev = _CELL_PREV_X_POS;
			} break;
			case _CELL_STEP_Z_POS: {
				next_z++;
				prev = _CELL_PREV_Z_NEG;
			} break;
			case _CELL_STEP_Z_NEG: {
				next_z--;
				prev = _CELL_PREV_Z_POS;
			} break;
			default:
				ERR_FAIL();
		}

		if (next_x < 0 || next_x >= len_x) {
			continue;
		}
		if (next_y < 0 || next_y >= len_y) {
			continue;
		}
		if (next_z < 0 || next_z >= len_z) {
			continue;
		}

		if (p_cell_status[next_x][next_y][next_z] & 3) {
			continue;
		}

		x = next_x;
		y = next_y;
		z = next_z;
		p_cell_status[x][y][z] |= prev;
	}
}

static inline void _build_faces(uint8_t ***p_cell_status, int x, int y, int z, int len_x, int len_y, int len_z, Vector<Face3> &p_faces) {
	ERR_FAIL_INDEX(x, len_x);
	ERR_FAIL_INDEX(y, len_y);
	ERR_FAIL_INDEX(z, len_z);

	if (p_cell_status[x][y][z] & _CELL_EXTERIOR) {
		return;
	}

#define vert(m_idx) Vector3(((m_idx)&4) >> 2, ((m_idx)&2) >> 1, (m_idx)&1)

	static const uint8_t indices[6][4] = {
		{ 7, 6, 4, 5 },
		{ 7, 3, 2, 6 },
		{ 7, 5, 1, 3 },
		{ 0, 2, 3, 1 },
		{ 0, 1, 5, 4 },
		{ 0, 4, 6, 2 },

	};

	for (int i = 0; i < 6; i++) {
		Vector3 face_points[4];
		int disp_x = x + ((i % 3) == 0 ? ((i < 3) ? 1 : -1) : 0);
		int disp_y = y + (((i - 1) % 3) == 0 ? ((i < 3) ? 1 : -1) : 0);
		int disp_z = z + (((i - 2) % 3) == 0 ? ((i < 3) ? 1 : -1) : 0);

		bool plot = false;

		if (disp_x < 0 || disp_x >= len_x) {
			plot = true;
		}
		if (disp_y < 0 || disp_y >= len_y) {
			plot = true;
		}
		if (disp_z < 0 || disp_z >= len_z) {
			plot = true;
		}

		if (!plot && (p_cell_status[disp_x][disp_y][disp_z] & _CELL_EXTERIOR)) {
			plot = true;
		}

		if (!plot) {
			continue;
		}

		for (int j = 0; j < 4; j++) {
			face_points[j] = vert(indices[i][j]) + Vector3(x, y, z);
		}

		p_faces.push_back(
				Face3(
						face_points[0],
						face_points[1],
						face_points[2]));

		p_faces.push_back(
				Face3(
						face_points[2],
						face_points[3],
						face_points[0]));
	}
}

Vector<Face3> Geometry3D::wrap_geometry(Vector<Face3> p_array, real_t *p_error) {
	int face_count = p_array.size();
	const Face3 *faces = p_array.ptr();
	constexpr double min_size = 1.0;
	constexpr int max_length = 20;

	AABB global_aabb;

	for (int i = 0; i < face_count; i++) {
		if (i == 0) {
			global_aabb = faces[i].get_aabb();
		} else {
			global_aabb.merge_with(faces[i].get_aabb());
		}
	}

	global_aabb.grow_by(0.01f); // Avoid numerical error.

	// Determine amount of cells in grid axis.
	int div_x, div_y, div_z;

	if (global_aabb.size.x / min_size < max_length) {
		div_x = (int)(global_aabb.size.x / min_size) + 1;
	} else {
		div_x = max_length;
	}

	if (global_aabb.size.y / min_size < max_length) {
		div_y = (int)(global_aabb.size.y / min_size) + 1;
	} else {
		div_y = max_length;
	}

	if (global_aabb.size.z / min_size < max_length) {
		div_z = (int)(global_aabb.size.z / min_size) + 1;
	} else {
		div_z = max_length;
	}

	Vector3 voxelsize = global_aabb.size;
	voxelsize.x /= div_x;
	voxelsize.y /= div_y;
	voxelsize.z /= div_z;

	// Create and initialize cells to zero.

	uint8_t ***cell_status = memnew_arr(uint8_t **, div_x);
	for (int i = 0; i < div_x; i++) {
		cell_status[i] = memnew_arr(uint8_t *, div_y);

		for (int j = 0; j < div_y; j++) {
			cell_status[i][j] = memnew_arr(uint8_t, div_z);

			for (int k = 0; k < div_z; k++) {
				cell_status[i][j][k] = 0;
			}
		}
	}

	// Plot faces into cells.

	for (int i = 0; i < face_count; i++) {
		Face3 f = faces[i];
		for (int j = 0; j < 3; j++) {
			f.vertex[j] -= global_aabb.position;
		}
		_plot_face(cell_status, 0, 0, 0, div_x, div_y, div_z, voxelsize, f);
	}

	// Determine which cells connect to the outside by traversing the outside and recursively flood-fill marking.

	for (int i = 0; i < div_x; i++) {
		for (int j = 0; j < div_y; j++) {
			_mark_outside(cell_status, i, j, 0, div_x, div_y, div_z);
			_mark_outside(cell_status, i, j, div_z - 1, div_x, div_y, div_z);
		}
	}

	for (int i = 0; i < div_z; i++) {
		for (int j = 0; j < div_y; j++) {
			_mark_outside(cell_status, 0, j, i, div_x, div_y, div_z);
			_mark_outside(cell_status, div_x - 1, j, i, div_x, div_y, div_z);
		}
	}

	for (int i = 0; i < div_x; i++) {
		for (int j = 0; j < div_z; j++) {
			_mark_outside(cell_status, i, 0, j, div_x, div_y, div_z);
			_mark_outside(cell_status, i, div_y - 1, j, div_x, div_y, div_z);
		}
	}

	// Build faces for the inside-outside cell divisors.

	Vector<Face3> wrapped_faces;

	for (int i = 0; i < div_x; i++) {
		for (int j = 0; j < div_y; j++) {
			for (int k = 0; k < div_z; k++) {
				_build_faces(cell_status, i, j, k, div_x, div_y, div_z, wrapped_faces);
			}
		}
	}

	// Transform face vertices to global coords.

	int wrapped_faces_count = wrapped_faces.size();
	Face3 *wrapped_faces_ptr = wrapped_faces.ptrw();

	for (int i = 0; i < wrapped_faces_count; i++) {
		for (int j = 0; j < 3; j++) {
			Vector3 &v = wrapped_faces_ptr[i].vertex[j];
			v = v * voxelsize;
			v += global_aabb.position;
		}
	}

	// clean up grid

	for (int i = 0; i < div_x; i++) {
		for (int j = 0; j < div_y; j++) {
			memdelete_arr(cell_status[i][j]);
		}

		memdelete_arr(cell_status[i]);
	}

	memdelete_arr(cell_status);
	if (p_error) {
		*p_error = voxelsize.length();
	}

	return wrapped_faces;
}

Geometry3D::MeshData Geometry3D::build_convex_mesh(const Vector<Plane> &p_planes) {
	MeshData mesh;

#define SUBPLANE_SIZE 1024.0

	real_t subplane_size = 1024.0; // Should compute this from the actual plane.
	for (int i = 0; i < p_planes.size(); i++) {
		Plane p = p_planes[i];

		Vector3 ref = Vector3(0.0, 1.0, 0.0);

		if (ABS(p.normal.dot(ref)) > 0.95f) {
			ref = Vector3(0.0, 0.0, 1.0); // Change axis.
		}

		Vector3 right = p.normal.cross(ref).normalized();
		Vector3 up = p.normal.cross(right).normalized();

		Vector3 center = p.center();

		// make a quad clockwise
		LocalVector<Vector3> vertices = {
			center - up * subplane_size + right * subplane_size,
			center - up * subplane_size - right * subplane_size,
			center + up * subplane_size - right * subplane_size,
			center + up * subplane_size + right * subplane_size
		};

		for (int j = 0; j < p_planes.size(); j++) {
			if (j == i) {
				continue;
			}

			LocalVector<Vector3> new_vertices;
			Plane clip = p_planes[j];

			if (clip.normal.dot(p.normal) > 0.95f) {
				continue;
			}

			if (vertices.size() < 3) {
				break;
			}

			for (uint32_t k = 0; k < vertices.size(); k++) {
				int k_n = (k + 1) % vertices.size();

				Vector3 edge0_A = vertices[k];
				Vector3 edge1_A = vertices[k_n];

				real_t dist0 = clip.distance_to(edge0_A);
				real_t dist1 = clip.distance_to(edge1_A);

				if (dist0 <= 0) { // Behind plane.

					new_vertices.push_back(vertices[k]);
				}

				// Check for different sides and non coplanar.
				if ((dist0 * dist1) < 0) {
					// Calculate intersection.
					Vector3 rel = edge1_A - edge0_A;

					real_t den = clip.normal.dot(rel);
					if (Math::is_zero_approx(den)) {
						continue; // Point too short.
					}

					real_t dist = -(clip.normal.dot(edge0_A) - clip.d) / den;
					Vector3 inters = edge0_A + rel * dist;
					new_vertices.push_back(inters);
				}
			}

			vertices = new_vertices;
		}

		if (vertices.size() < 3) {
			continue;
		}

		// Result is a clockwise face.

		MeshData::Face face;

		// Add face indices.
		for (uint32_t j = 0; j < vertices.size(); j++) {
			int idx = -1;
			for (uint32_t k = 0; k < mesh.vertices.size(); k++) {
				if (mesh.vertices[k].distance_to(vertices[j]) < 0.001f) {
					idx = k;
					break;
				}
			}

			if (idx == -1) {
				idx = mesh.vertices.size();
				mesh.vertices.push_back(vertices[j]);
			}

			face.indices.push_back(idx);
		}
		face.plane = p;
		mesh.faces.push_back(face);

		// Add edge.

		for (uint32_t j = 0; j < face.indices.size(); j++) {
			int a = face.indices[j];
			int b = face.indices[(j + 1) % face.indices.size()];

			bool found = false;
			int found_idx = -1;
			for (uint32_t k = 0; k < mesh.edges.size(); k++) {
				if (mesh.edges[k].vertex_a == a && mesh.edges[k].vertex_b == b) {
					found = true;
					found_idx = k;
					break;
				}
				if (mesh.edges[k].vertex_b == a && mesh.edges[k].vertex_a == b) {
					found = true;
					found_idx = k;
					break;
				}
			}

			if (found) {
				mesh.edges[found_idx].face_b = j;
				continue;
			}
			MeshData::Edge edge;
			edge.vertex_a = a;
			edge.vertex_b = b;
			edge.face_a = j;
			edge.face_b = -1;
			mesh.edges.push_back(edge);
		}
	}

	return mesh;
}

Vector<Plane> Geometry3D::build_box_planes(const Vector3 &p_extents) {
	Vector<Plane> planes = {
		Plane(Vector3(1, 0, 0), p_extents.x),
		Plane(Vector3(-1, 0, 0), p_extents.x),
		Plane(Vector3(0, 1, 0), p_extents.y),
		Plane(Vector3(0, -1, 0), p_extents.y),
		Plane(Vector3(0, 0, 1), p_extents.z),
		Plane(Vector3(0, 0, -1), p_extents.z)
	};

	return planes;
}

Vector<Plane> Geometry3D::build_cylinder_planes(real_t p_radius, real_t p_height, int p_sides, Vector3::Axis p_axis) {
	ERR_FAIL_INDEX_V(p_axis, 3, Vector<Plane>());

	Vector<Plane> planes;

	const double sides_step = Math_TAU / p_sides;
	for (int i = 0; i < p_sides; i++) {
		Vector3 normal;
		normal[(p_axis + 1) % 3] = Math::cos(i * sides_step);
		normal[(p_axis + 2) % 3] = Math::sin(i * sides_step);

		planes.push_back(Plane(normal, p_radius));
	}

	Vector3 axis;
	axis[p_axis] = 1.0;

	planes.push_back(Plane(axis, p_height * 0.5f));
	planes.push_back(Plane(-axis, p_height * 0.5f));

	return planes;
}

Vector<Plane> Geometry3D::build_sphere_planes(real_t p_radius, int p_lats, int p_lons, Vector3::Axis p_axis) {
	ERR_FAIL_INDEX_V(p_axis, 3, Vector<Plane>());

	Vector<Plane> planes;

	Vector3 axis;
	axis[p_axis] = 1.0;

	Vector3 axis_neg;
	axis_neg[(p_axis + 1) % 3] = 1.0;
	axis_neg[(p_axis + 2) % 3] = 1.0;
	axis_neg[p_axis] = -1.0;

	const double lon_step = Math_TAU / p_lons;
	for (int i = 0; i < p_lons; i++) {
		Vector3 normal;
		normal[(p_axis + 1) % 3] = Math::cos(i * lon_step);
		normal[(p_axis + 2) % 3] = Math::sin(i * lon_step);

		planes.push_back(Plane(normal, p_radius));

		for (int j = 1; j <= p_lats; j++) {
			Vector3 plane_normal = normal.lerp(axis, j / (real_t)p_lats).normalized();
			planes.push_back(Plane(plane_normal, p_radius));
			planes.push_back(Plane(plane_normal * axis_neg, p_radius));
		}
	}

	return planes;
}

Vector<Plane> Geometry3D::build_capsule_planes(real_t p_radius, real_t p_height, int p_sides, int p_lats, Vector3::Axis p_axis) {
	ERR_FAIL_INDEX_V(p_axis, 3, Vector<Plane>());

	Vector<Plane> planes;

	Vector3 axis;
	axis[p_axis] = 1.0;

	Vector3 axis_neg;
	axis_neg[(p_axis + 1) % 3] = 1.0;
	axis_neg[(p_axis + 2) % 3] = 1.0;
	axis_neg[p_axis] = -1.0;

	const double sides_step = Math_TAU / p_sides;
	for (int i = 0; i < p_sides; i++) {
		Vector3 normal;
		normal[(p_axis + 1) % 3] = Math::cos(i * sides_step);
		normal[(p_axis + 2) % 3] = Math::sin(i * sides_step);

		planes.push_back(Plane(normal, p_radius));

		for (int j = 1; j <= p_lats; j++) {
			Vector3 plane_normal = normal.lerp(axis, j / (real_t)p_lats).normalized();
			Vector3 position = axis * p_height * 0.5f + plane_normal * p_radius;
			planes.push_back(Plane(plane_normal, position));
			planes.push_back(Plane(plane_normal * axis_neg, position * axis_neg));
		}
	}

	return planes;
}

Vector<Vector3> Geometry3D::compute_convex_mesh_points(const Plane *p_planes, int p_plane_count) {
	Vector<Vector3> points;

	// Iterate through every unique combination of any three planes.
	for (int i = p_plane_count - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			for (int k = j - 1; k >= 0; k--) {
				// Find the point where these planes all cross over (if they
				// do at all).
				Vector3 convex_shape_point;
				if (p_planes[i].intersect_3(p_planes[j], p_planes[k], &convex_shape_point)) {
					// See if any *other* plane excludes this point because it's
					// on the wrong side.
					bool excluded = false;
					for (int n = 0; n < p_plane_count; n++) {
						if (n != i && n != j && n != k) {
							real_t dp = p_planes[n].normal.dot(convex_shape_point);
							if (dp - p_planes[n].d > (real_t)CMP_EPSILON) {
								excluded = true;
								break;
							}
						}
					}

					// Only add the point if it passed all tests.
					if (!excluded) {
						points.push_back(convex_shape_point);
					}
				}
			}
		}
	}

	return points;
}

#define square(m_s) ((m_s) * (m_s))
#define INF 1e20

/* dt of 1d function using squared distance */
static void edt(float *f, int stride, int n) {
	float *d = (float *)alloca(sizeof(float) * n + sizeof(int) * n + sizeof(float) * (n + 1));
	int *v = reinterpret_cast<int *>(&(d[n]));
	float *z = reinterpret_cast<float *>(&v[n]);

	int k = 0;
	v[0] = 0;
	z[0] = -INF;
	z[1] = +INF;
	for (int q = 1; q <= n - 1; q++) {
		float s = ((f[q * stride] + square(q)) - (f[v[k] * stride] + square(v[k]))) / (2 * q - 2 * v[k]);
		while (s <= z[k]) {
			k--;
			s = ((f[q * stride] + square(q)) - (f[v[k] * stride] + square(v[k]))) / (2 * q - 2 * v[k]);
		}
		k++;
		v[k] = q;

		z[k] = s;
		z[k + 1] = +INF;
	}

	k = 0;
	for (int q = 0; q <= n - 1; q++) {
		while (z[k + 1] < q) {
			k++;
		}
		d[q] = square(q - v[k]) + f[v[k] * stride];
	}

	for (int i = 0; i < n; i++) {
		f[i * stride] = d[i];
	}
}

#undef square

Vector<uint32_t> Geometry3D::generate_edf(const Vector<bool> &p_voxels, const Vector3i &p_size, bool p_negative) {
	uint32_t float_count = p_size.x * p_size.y * p_size.z;

	ERR_FAIL_COND_V((uint32_t)p_voxels.size() != float_count, Vector<uint32_t>());

	float *work_memory = memnew_arr(float, float_count);
	for (uint32_t i = 0; i < float_count; i++) {
		work_memory[i] = INF;
	}

	uint32_t y_mult = p_size.x;
	uint32_t z_mult = y_mult * p_size.y;

	//plot solid cells
	{
		const bool *voxr = p_voxels.ptr();
		for (uint32_t i = 0; i < float_count; i++) {
			bool plot = voxr[i];
			if (p_negative) {
				plot = !plot;
			}
			if (plot) {
				work_memory[i] = 0;
			}
		}
	}

	//process in each direction

	//xy->z

	for (int i = 0; i < p_size.x; i++) {
		for (int j = 0; j < p_size.y; j++) {
			edt(&work_memory[i + j * y_mult], z_mult, p_size.z);
		}
	}

	//xz->y

	for (int i = 0; i < p_size.x; i++) {
		for (int j = 0; j < p_size.z; j++) {
			edt(&work_memory[i + j * z_mult], y_mult, p_size.y);
		}
	}

	//yz->x
	for (int i = 0; i < p_size.y; i++) {
		for (int j = 0; j < p_size.z; j++) {
			edt(&work_memory[i * y_mult + j * z_mult], 1, p_size.x);
		}
	}

	Vector<uint32_t> ret;
	ret.resize(float_count);
	{
		uint32_t *w = ret.ptrw();
		for (uint32_t i = 0; i < float_count; i++) {
			w[i] = uint32_t(Math::sqrt(work_memory[i]));
		}
	}

	memdelete_arr(work_memory);

	return ret;
}

Vector<int8_t> Geometry3D::generate_sdf8(const Vector<uint32_t> &p_positive, const Vector<uint32_t> &p_negative) {
	ERR_FAIL_COND_V(p_positive.size() != p_negative.size(), Vector<int8_t>());
	Vector<int8_t> sdf8;
	int s = p_positive.size();
	sdf8.resize(s);

	const uint32_t *rpos = p_positive.ptr();
	const uint32_t *rneg = p_negative.ptr();
	int8_t *wsdf = sdf8.ptrw();
	for (int i = 0; i < s; i++) {
		int32_t diff = int32_t(rpos[i]) - int32_t(rneg[i]);
		wsdf[i] = CLAMP(diff, -128, 127);
	}
	return sdf8;
}

/*
 * Line-Circle distance (based on DistLine3Circle3.h from David Eberly's Geometric Tools)
 */

struct LineCirclePointPair {
	Vector3 line_point;
	Vector3 circle_point;
	real_t squared_distance = 0.0;
	bool equidistant = false;
};

struct LineCirclePointPairSort {
	bool operator()(const LineCirclePointPair p_a, const LineCirclePointPair p_b) const {
		return p_a.squared_distance < p_b.squared_distance;
	}
};

void get_point_pair_from_line_parameter(const Vector3 &p_circle_center, const Vector3 &p_circle_normal, const real_t p_circle_radius, const Vector3 &p_line_direction, const Vector3 &p_D, const real_t p_t, Vector3 &r_line_point, Vector3 &r_circle_point) {
	// GetPair by David Eberly
	Vector3 delta = p_D + p_t * p_line_direction;
	r_line_point = p_circle_center + delta;
	delta -= p_circle_normal.dot(delta) * p_circle_normal;
	delta.normalize();
	r_circle_point = p_circle_center + p_circle_radius * delta;
}

Vector3 get_some_orthogonal_vector(const Vector3 &p_vector, bool p_unit_length) {
	// GetOrthogonal by David Eberly

	// Construct a single vector orthogonal to the nonzero input vector. If
	// the maximum absolute component occurs at index i, then the orthogonal
	// vector U has u[i] = v[i+1], u[i+1] = -v[i], and all other components
	// zero. The index addition i+1 is computed modulo N. If the input vector
	// is zero, the output vector is zero. If the input vector is empty, the
	// output vector is empty. If you want the output to be a unit-length
	// vector, set unitLength to 'true'.

	real_t cmax = 0.0;
	size_t imax = 0;
	for (size_t i = 0; i < 3; ++i) {
		real_t c = Math::abs(p_vector[i]);
		if (c > cmax) {
			cmax = c;
			imax = i;
		}
	}

	Vector3 result;
	if (cmax > 0.0) {
		size_t inext = imax + 1;
		if (inext == 3) {
			inext = 0;
		}
		result[imax] = p_vector[inext];
		result[inext] = -p_vector[imax];
		if (p_unit_length) {
			result.normalize();
		}
	}
	return result;
}

void Geometry3D::get_closest_points_between_line_and_circle(const Vector3 &p_line_origin, const Vector3 &p_line_direction, const Vector3 &p_circle_center, const Vector3 &p_circle_normal, const real_t p_circle_radius, Vector<Vector3> &r_line_closest, Vector<Vector3> &r_circle_closest, size_t &r_num_closest_pairs, bool &r_equidistant) {
	// Based on David Eberly's Distance Between a Line and a Circle: polynomial-based algorithm.
	// The only difference is the use of Godot's types, and the use of a numerical iterative root finder.

	Vector3 *circle_closest = r_circle_closest.ptrw();
	Vector3 *line_closest = r_line_closest.ptrw();

	Vector3 D = p_line_origin - p_circle_center;
	Vector3 NxM = p_circle_normal.cross(p_line_direction);
	Vector3 NxD = p_circle_normal.cross(D);
	real_t t;

	if (NxM != Vector3(0, 0, 0)) {
		if (NxD != Vector3(0, 0, 0)) {
			real_t NdM = p_circle_normal.dot(p_line_direction);
			if (NdM != 0.0) {
				// H(t) = (a*t^2 + 2*b*t + c)*(t + d)^2
				//        - r^2*(a*t + b)^2
				//      = h0 + h1*t + h2*t^2 + h3*t^3 + h4*t^4
				real_t a = NxM.dot(NxM);
				real_t b = NxM.dot(NxD);
				real_t c = NxD.dot(NxD);
				real_t d = p_line_direction.dot(D);
				real_t rSqr = p_circle_radius * p_circle_radius;
				real_t aSqr = a * a;
				real_t bSqr = b * b;
				real_t dSqr = d * d;
				real_t h0 = c * dSqr - bSqr * rSqr;
				real_t h1 = 2.0 * (c * d + b * dSqr - a * b * rSqr);
				real_t h2 = c + 4.0 * b * d + a * dSqr - aSqr * rSqr;
				real_t h3 = 2.0 * (b + a * d);
				real_t h4 = a;

				Vector<real_t> real_roots = real_roots_of_nonconstant_polynomial({ h4, h3, h2, h1, h0 });

				Vector<LineCirclePointPair> point_pairs;
				for (real_t real_root : real_roots) {
					t = real_root;
					LineCirclePointPair point_pair;
					Vector3 NxDelta = NxD + t * NxM;
					if (NxDelta != Vector3(0, 0, 0)) {
						get_point_pair_from_line_parameter(p_circle_center, p_circle_normal, p_circle_radius, p_line_direction, D, t, point_pair.line_point, point_pair.circle_point);
						point_pair.equidistant = false;
					} else {
						Vector3 U = get_some_orthogonal_vector(p_circle_normal, true);
						point_pair.line_point = p_circle_center;
						point_pair.circle_point = p_circle_center + p_circle_radius * U;
						point_pair.equidistant = true;
					}
					Vector3 diff = point_pair.line_point - point_pair.circle_point;
					point_pair.squared_distance = diff.dot(diff);
					point_pairs.push_back(point_pair);
				}

				if (point_pairs.size() == 0) { // Failed to find any real roots.
					// TODO: investigate this.
					return;
				}

				point_pairs.sort_custom<LineCirclePointPairSort>();

				r_num_closest_pairs = 1;
				circle_closest[0] = point_pairs[0].circle_point;
				line_closest[0] = point_pairs[0].line_point;
				if (point_pairs.size() > 1 && Math::is_equal_approx(point_pairs[0].squared_distance, point_pairs[1].squared_distance)) {
					r_num_closest_pairs = 2;
					circle_closest[1] = point_pairs[1].circle_point;
					line_closest[1] = point_pairs[1].line_point;
				}
			} else {
				// The line is parallel to the plane of the circle.
				// The polynomial has the form
				// H(t) = (t+v)^2*[(t+v)^2-(r^2-u^2)].
				real_t u = NxM.dot(D);
				real_t v = p_line_direction.dot(D);
				real_t discr = p_circle_radius * p_circle_radius - u * u;
				if (discr > 0.0) {
					r_num_closest_pairs = 2;
					real_t rootDiscr = Math::sqrt(discr);
					t = -v + rootDiscr;
					get_point_pair_from_line_parameter(p_circle_center, p_circle_normal, p_circle_radius, p_line_direction, D, t, line_closest[0], circle_closest[0]);
					t = -v - rootDiscr;
					get_point_pair_from_line_parameter(p_circle_center, p_circle_normal, p_circle_radius, p_line_direction, D, t, line_closest[1], circle_closest[1]);
				} else {
					r_num_closest_pairs = 1;
					t = -v;
					get_point_pair_from_line_parameter(p_circle_center, p_circle_normal, p_circle_radius, p_line_direction, D, t, line_closest[0], circle_closest[0]);
				}
			}
		} else {
			// The line is C+t*M, where M is not parallel to N. The
			// polynomial is
			// H(t) = |Cross(N,M)|^2*t^2*(t^2 - r^2*|Cross(N,M)|^2)
			// where root t = 0 does not correspond to the global
			// minimum. The other roots produce the global minimum.
			r_num_closest_pairs = 2;
			t = p_circle_radius * NxM.length();
			get_point_pair_from_line_parameter(p_circle_center, p_circle_normal, p_circle_radius, p_line_direction, D, t, line_closest[0], circle_closest[0]);
			t = -t;
			get_point_pair_from_line_parameter(p_circle_center, p_circle_normal, p_circle_radius, p_line_direction, D, t, line_closest[1], circle_closest[1]);
		}
		r_equidistant = false;
	} else {
		if (NxD != Vector3(0, 0, 0)) {
			// The line is A+t*N (perpendicular to plane) but with
			// A != C. The polyhomial is
			// H(t) = |Cross(N,D)|^2*(t + Dot(M,D))^2.
			r_num_closest_pairs = 1;
			t = -p_line_direction.dot(D);
			get_point_pair_from_line_parameter(p_circle_center, p_circle_normal, p_circle_radius, p_line_direction, D, t, line_closest[0], circle_closest[0]);
			r_equidistant = false;
		} else {
			// The line is C+t*N, so C is the closest point for the
			// line and all circle points are equidistant from it.
			Vector3 U = get_some_orthogonal_vector(p_circle_normal, true);
			r_num_closest_pairs = 1;
			line_closest[0] = p_circle_center;
			circle_closest[0] = p_circle_center + p_circle_radius * U;
			r_equidistant = true;
		}
	}
}

/*
 * Circle-Circle distance (based on DistCircle3Circle3.h from David Eberly's Geometric Tools)
 */

struct CircleCirclePointPair {
	Vector3 circle0_point;
	Vector3 circle1_point;
	real_t squared_distance = 0.0;
	bool equidistant = false;
};

struct CircleCirclePointPairSort {
	bool operator()(const CircleCirclePointPair p_a, const CircleCirclePointPair p_b) const {
		return p_a.squared_distance < p_b.squared_distance;
	}
};

bool compute_orthonormal_basis(Vector3 &r_vector0, Vector3 &r_vector1, Vector3 &r_vector2) {
	// ComputeOrthonormalBasis by David Eberly.
	// Compute a right-handed orthonormal basis from one nonzero vector.
	// The function returns true when the basis is computed successfully,
	// in which case the matrix [v0 v1 v2] is a rotation matrix.
	// If the function returns false, the outputs v0, v1 and v2 are invalid.

	r_vector0.normalize();
	if (r_vector0 == Vector3(0, 0, 0)) {
		r_vector1 = Vector3(0, 0, 0);
		r_vector2 = Vector3(0, 0, 0);
		return false;
	}
	if (Math::abs(r_vector0[0]) > Math::abs(r_vector0[1])) {
		r_vector1 = Vector3(-r_vector0[2], 0.0, r_vector0[0]);
	} else {
		r_vector1 = Vector3(0.0, r_vector0[2], -r_vector0[1]);
	}
	r_vector1.normalize();
	if (r_vector1 == Vector3(0, 0, 0)) {
		r_vector2 = Vector3(0, 0, 0);
		return false;
	}
	r_vector2 = r_vector0.cross(r_vector1);
	return r_vector2 != Vector3(0, 0, 0);
}

void Geometry3D::get_closest_points_between_circle_and_circle(const Vector3 &p_circle0_center, const Vector3 &p_circle0_normal, const real_t p_circle0_radius, const Vector3 &p_circle1_center, const Vector3 &p_circle1_normal, const real_t p_circle1_radius, Vector<Vector3> &r_circle0_closest, Vector<Vector3> &r_circle1_closest, size_t &r_num_closest_pairs, bool &r_equidistant) {
	// Based on David Eberly's Distance Between Two Circles algorithm.
	// The only difference is the use of Godot's types (and a minimalist Polynomial type), and the use of a different numerical iterative root finder.

	Vector3 *circle0_closest = r_circle0_closest.ptrw();
	Vector3 *circle1_closest = r_circle1_closest.ptrw();

	Vector3 N0 = p_circle0_normal;
	Vector3 N1 = p_circle1_normal;
	real_t r0 = p_circle0_radius;
	real_t r1 = p_circle1_radius;
	Vector3 D = p_circle1_center - p_circle0_center;
	Vector3 N0xN1 = N0.cross(N1);

	if (N0xN1 != Vector3(0, 0, 0)) // comparison to the zero vector
	{
		// Get parameters for constructing the degree-8 polynomial phi.
		real_t r0sqr = r0 * r0;
		real_t r1sqr = r1 * r1;

		// Compute U1 and V1 for the plane of circle1.
		Vector3 U1, V1;
		compute_orthonormal_basis(N1, U1, V1);

		// Construct the polynomial phi(cos(theta)).
		Vector3 N0xD = N0.cross(D);
		Vector3 N0xU1 = N0.cross(U1);
		Vector3 N0xV1 = N0.cross(V1);
		real_t a0 = r1 * D.dot(U1);
		real_t a1 = r1 * D.dot(V1);
		real_t a2 = N0xD.dot(N0xD);
		real_t a3 = r1 * N0xD.dot(N0xU1);
		real_t a4 = r1 * N0xD.dot(N0xV1);
		real_t a5 = r1sqr * N0xU1.dot(N0xU1);
		real_t a6 = r1sqr * N0xU1.dot(N0xV1);
		real_t a7 = r1sqr * N0xV1.dot(N0xV1);
		Polynomial p0({ a2 + a7, real_t(2.0) * a3, a5 - a7 });
		Polynomial p1({ real_t(2.0) * a4, real_t(2.0) * a6 });
		Polynomial p2({ real_t(0.0), a1 });
		Polynomial p3(Vector<real_t>{ -a0 });
		Polynomial p4({ -a6, a4, real_t(2.0) * a6 });
		Polynomial p5({ -a3, a7 - a5 });
		Polynomial tmp0({ real_t(1.0), real_t(0.0), real_t(-1.0) });
		Polynomial tmp1 = p2 * p2 + tmp0 * p3 * p3;
		Polynomial tmp2 = real_t(2.0) * p2 * p3;
		Polynomial tmp3 = p4 * p4 + tmp0 * p5 * p5;
		Polynomial tmp4 = real_t(2.0) * p4 * p5;
		Polynomial p6 = p0 * tmp1 + tmp0 * p1 * tmp2 - r0sqr * tmp3;
		Polynomial p7 = p0 * tmp2 + p1 * tmp1 - r0sqr * tmp4;

		std::array<std::pair<real_t, real_t>, 16> pairs{}; // TODO: no std
		size_t numPairs = 0;
		if (p7.degree() > 0 || p7[0] != 0.0) {
			// H(cs,sn) = p6(cs) + sn * p7(cs)
			Polynomial phi = p6 * p6 - tmp0 * p7 * p7;

			ERR_FAIL_COND_MSG(phi.degree() == 0, "Unexpected degree for phi.");

			Vector<real_t> real_roots = real_roots_of_nonconstant_polynomial(phi.coefficients_decreasing());
			// TODO: unique roots?

			for (real_t const &cs : real_roots) {
				if (Math::abs(cs) <= 1.0) {
					real_t temp = p7(cs);
					if (temp != 0.0) {
						real_t sn = -p6(cs) / temp;
						pairs[numPairs++] = std::make_pair(cs, sn);
					} else {
						temp = MAX(1.0 - cs * cs, 0.0);
						real_t sn = Math::sqrt(temp);
						pairs[numPairs++] = std::make_pair(cs, sn);
						if (sn != 0.0) {
							pairs[numPairs++] = std::make_pair(cs, -sn);
						}
					}
				}
			}
		} else {
			// H(cs,sn) = p6(cs)
			ERR_FAIL_COND_MSG(p6.degree() == 0, "Unexpected degree for p6.");

			Vector<real_t> real_roots = real_roots_of_nonconstant_polynomial(p6.coefficients_decreasing());
			// TODO: unique roots?

			for (real_t const &cs : real_roots) {
				if (Math::abs(cs) <= 1.0) {
					real_t temp = MAX(1.0 - cs * cs, 0.0);
					real_t sn = Math::sqrt(temp);
					pairs[numPairs++] = std::make_pair(cs, sn);
					if (sn != 0.0) {
						pairs[numPairs++] = std::make_pair(cs, -sn);
					}
				}
			}
		}

		if (numPairs == 0) { // Failed to find any real roots.
			return; // TODO: investigate this.
		}

		Vector<CircleCirclePointPair> candidates;
		candidates.resize_zeroed(numPairs); // TODO: correct?
		for (size_t i = 0; i < numPairs; ++i) {
			CircleCirclePointPair &info = candidates.ptrw()[i];
			Vector3 delta = D + r1 * (pairs[i].first * U1 + pairs[i].second * V1);
			info.circle1_point = p_circle0_center + delta;
			real_t N0dDelta = N0.dot(delta);
			real_t lenN0xDelta = N0.cross(delta).length();
			if (lenN0xDelta > 0.0) {
				real_t diff = lenN0xDelta - r0;
				info.squared_distance = N0dDelta * N0dDelta + diff * diff;
				delta -= N0dDelta * p_circle0_normal;
				delta.normalize();
				info.circle0_point = p_circle0_center + r0 * delta;
				info.equidistant = false;
			} else {
				Vector3 r0U0 = r0 * get_some_orthogonal_vector(N0, true);
				Vector3 diff = delta - r0U0;
				info.squared_distance = diff.dot(diff);
				info.circle0_point = p_circle0_center + r0U0;
				info.equidistant = true;
			}
		}

		candidates.sort_custom<CircleCirclePointPairSort>();

		r_num_closest_pairs = 1;
		circle0_closest[0] = candidates[0].circle0_point;
		circle1_closest[0] = candidates[0].circle1_point;
		r_equidistant = candidates[0].equidistant;
		if (candidates.size() > 1 && candidates[1].squared_distance == candidates[0].squared_distance) {
			r_num_closest_pairs = 2;
			circle0_closest[1] = candidates[1].circle0_point;
			circle1_closest[1] = candidates[1].circle1_point;
		}
	} else {
		// The planes of the circles are parallel. Whether the planes
		// are the same or different, the problem reduces to
		// determining how two circles in the same plane are
		// separated, tangent with one circle outside the other,
		// overlapping or one circle contained inside the other
		// circle.

		// Based on David Eberly's DoQueryParallelPlanes.

		real_t N0dD = p_circle0_normal.dot(D);
		Vector3 normProj = N0dD * p_circle0_normal;
		Vector3 compProj = D - normProj;
		Vector3 U = compProj;
		real_t d = U.length();
		U.normalize();

		// The configuration is determined by the relative location of the
		// intervals of projection of the circles on to the D-line.
		// Circle0 projects to [-r0,r0] and circle1 projects to
		// [d-r1,d+r1].
		real_t r0 = p_circle0_radius;
		real_t r1 = p_circle1_radius;
		real_t dmr1 = d - r1;
		//real_t distance;
		if (dmr1 >= r0) // d >= r0 + r1
		{
			// The circles are separated (d > r0 + r1) or tangent with one
			// outside the other (d = r0 + r1).
			//distance = dmr1 - r0;
			r_num_closest_pairs = 1;
			circle0_closest[0] = p_circle0_center + r0 * U;
			circle1_closest[0] = p_circle1_center - r1 * U;
			r_equidistant = false;
		} else { // d < r0 + r1
			// The cases implicitly use the knowledge that d >= 0.
			real_t dpr1 = d + r1;
			if (dpr1 <= r0) {
				// Circle1 is inside circle0.
				//distance = r0 - dpr1;
				r_num_closest_pairs = 1;
				if (d > 0.0) {
					circle0_closest[0] = p_circle0_center + r0 * U;
					circle1_closest[0] = p_circle1_center + r1 * U;
					r_equidistant = false;
				} else {
					// The circles are concentric, so U = (0,0,0).
					// Construct a vector perpendicular to N0 to use for
					// closest points.
					U = get_some_orthogonal_vector(p_circle0_normal, true);
					circle0_closest[0] = p_circle0_center + r0 * U;
					circle1_closest[0] = p_circle1_center + r1 * U;
					r_equidistant = true;
				}
			} else if (dmr1 <= -r0) {
				// Circle0 is inside circle1.
				//distance = -r0 - dmr1;
				r_num_closest_pairs = 1;
				if (d > 0.0) {
					circle0_closest[0] = p_circle0_center - r0 * U;
					circle1_closest[0] = p_circle1_center - r1 * U;
					r_equidistant = false;
				} else {
					// The circles are concentric, so U = (0,0,0).
					// Construct a vector perpendicular to N0 to use for
					// closest points.
					U = get_some_orthogonal_vector(p_circle0_normal, true);
					circle0_closest[0] = p_circle0_center + r0 * U;
					circle1_closest[0] = p_circle1_center + r1 * U;
					r_equidistant = true;
				}
			} else {
				// The circles are overlapping.  The two points of
				// intersection are C0 + s*(C1-C0) +/- h*Cross(N,U), where
				// s = (1 + (r0^2 - r1^2)/d^2)/2 and
				// h = sqrt(r0^2 - s^2 * d^2).
				real_t r0sqr = r0 * r0;
				real_t r1sqr = r1 * r1;
				real_t dsqr = d * d;
				real_t s = (1.0 + (r0sqr - r1sqr) / dsqr) / 2.0;
				real_t arg = MAX(r0sqr - dsqr * s * s, 0.0);
				real_t h = Math::sqrt(arg);
				Vector3 midpoint = p_circle0_center + s * compProj;
				Vector3 hNxU = h * p_circle0_normal.cross(U);
				//distance = 0.0;
				r_num_closest_pairs = 2;
				circle0_closest[0] = midpoint + hNxU;
				circle0_closest[1] = midpoint - hNxU;
				circle1_closest[0] = circle0_closest[0] + normProj;
				circle1_closest[1] = circle0_closest[1] + normProj;
				r_equidistant = false;
			}
		}
	}
}
