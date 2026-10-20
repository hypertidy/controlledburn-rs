# Topology and meshes: the surrounding software

Status: discussion. Companion to `topology-mesh.md`. Not exhaustive: enough
to show where constrained triangulation, arrangements and shared-edge
models are reachable from, and what controlledburn tooling around them
could look like. Most of the engine comparison comes from the
trianglewins harness (github.com/mdsumner/trianglewins, October 2026),
which runs the same inputs through four constrained triangulators and
judges them with code that belongs to none of them.

## Why it matters here

`topology-mesh.md` proposes a shared-arc input (arcs with a left and a
right id) and an overlay output (boundary pieces per cell with left and
right ids). Both assume someone else has already built the topology:
noded the segments, merged shared vertices and edges, and labelled the
faces. The tools below are where that happens. controlledburn should
take their output, not reimplement them.

## The shared input: a planar straight line graph

The engines that matter all accept the same thing: vertices (x, y) and
segments as index pairs (s0, s1), a PSLG. Polygons, lines and points are
all PSLGs. The hypertidy packages (cdtr, laridae, trowel) share this
contract, plus three outputs that are directly useful to a burn:

- a vertex table with `origin` (input, crossing, steiner), so new
  vertices are traceable;
- a triangle table with **depth** (constraint crossings from outside),
  which decides inside and outside topologically and works for a
  coverage with shared edges, where ring parity and hole seeds do not;
- a segment table mapping each mesh edge back to the input segments it
  came from, with a `count` when a boundary was given more than once.

The triangle form is written up as the hypertidy mesh spec
(github.com/hypertidy/meshcore, mesh-spec v0.1); draft work in
trianglewins also covers a `geoarrow.triangles` extension type and a
UGRID-in-Zarr pilot read by the Python stack.

From a labelled triangle mesh, the arcs for controlledburn's shared-arc
mode are the constrained edges whose two triangles carry different
feature labels (or one triangle and the outside), with those labels as
left and right ids. Chaining them between nodes of degree other than two
gives arcs rather than single segments.

## Engines

| engine | kind | language | licence | notes |
|---|---|---|---|---|
| CGAL `Constrained_triangulation_plus_2`, `Delaunay_mesher_2` | CDT, refinement, editable | C++ (header only) | GPL for Mesh_2 and Triangulation_2 (commercial licence from GeometryFactory) | exact predicates (Epick); remembers input constraints after splits; insert and remove at any time |
| CGAL `Arrangement_2` | planar overlay with full face, edge and vertex incidence (DCEL) | C++ | GPL | the reference for "overlay with left and right faces"; not exercised by trianglewins |
| CDT (artem-ogre) | CDT, conforming, refinement | C++ (header only) | MPL-2.0 | permissive; refinement needs a closed boundary |
| Triangle (Shewchuk) | CDT, quality refinement | C | non-commercial | the reference Ruppert mesher, fastest refinement in the harness; refuses duplicate vertices |
| spade | CDT, refinement | Rust | MIT or Apache-2.0 | used by trowel; counts a doubly given edge once in depth |
| earcut (Mapbox) | ear cutting, one polygon at a time | C++, JS, Rust ports | ISC | for drawing; no shared edges, no Steiner points, poor angles |
| GEOS polygon triangulation | ear cutting plus flips, per polygon | C++ (C API) | LGPL | identical triangles to the CDT on valid coverages, but geometry out, no topology; lines dropped |
| Qhull | point Delaunay | C | Qhull licence | points only |

What the harness found, in short: on the shared contract every backend
kept every input segment and vertex, was constrained Delaunay, and gave
the same triangle count and depth where no boundary was given twice.
Depth differs when one is: laridae and cdtr count an edge covered by k
input segments k times, trowel once. A consumer has to name which.

## Where each is reachable from

| | R | Python | Rust |
|---|---|---|---|
| CGAL | laridae (via RcppCGAL headers, patched for R) | CGAL SWIG bindings (not checked here) | no maintained binding found; through C++ |
| CDT | cdtr | PythonCDT (by the CDT author) | through C++ |
| Triangle | RTriangle | triangle, meshpy | through C |
| spade | trowel (extendr) | | spade directly; geo's spade triangulation |
| earcut | decido | mapbox-earcut | earcutr; geo's earcut triangulation |
| GEOS | geos, sf, terra, gdalraster | shapely | geos crate |
| Qhull | geometry, tessellation | scipy.spatial | |

For a Rust crate the practical choices are spade (pure Rust, permissive,
already used through trowel) for triangulation, and C++ through a thin
FFI for CGAL or CDT where their guarantees (exact constraint hierarchy,
arrangements, editing) are needed. geo covers the simple-features side:
`Relate` for exact DE-9IM, boolean operations, and spade and earcut
triangulation behind traits.

## What the tooling could look like

Three layers, each usable alone:

1. **Topology builders** (outside controlledburn): any of the engines
   above, or a GEOS or geo noding step, turn features into a PSLG with
   shared vertices, then into labelled faces (depth or seeds).
2. **An arcs adapter** (small, could live in wxylib or an examples
   module): reads the shared tables (vertices, segments, triangles with
   labels) and emits arcs with left and right ids for the burn. It
   needs no geometry beyond chaining edges.
3. **controlledburn**: burns arcs once per edge, and optionally returns
   the overlay pieces. Downstream, the DE-9IM cell tables,
   ALL_TOUCHED masks and conservative regridding weights read those
   pieces.

Going the other way, the overlay pieces plus cell runs are themselves a
PSLG (pieces as segments, cell walls as segments). Fed back to a
constrained triangulator, that gives a mesh of polygon-in-cell faces,
which is the mesh-to-grid overlay used for conservative remapping (ESMF
builds the same thing internally). That round trip is also a cheap
consistency test for the overlay output.

## Not covered

- Spherical and geodesic edges (CGAL has spherical triangulations; the
  edge model question is the same one raised for bigcurve).
- 3D and 2.5D (TINs, surface meshes).
- Unstructured ocean and climate model meshes beyond the UGRID pointer.
