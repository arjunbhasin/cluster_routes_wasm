use wasm_bindgen::prelude::*;
use petgraph::graph::UnGraph;
use petgraph::algo::kosaraju_scc;
use std::collections::HashMap;

#[wasm_bindgen]
pub fn cluster_routes_js(routes_data: JsValue, threshold: f64) -> Result<JsValue, JsValue>{
    // Deserialize routes_data from JsValue into Vec<Vec<[f64; 2]>>
    let routes_data: Vec<Vec<[f64; 2]>> = serde_wasm_bindgen::from_value(routes_data)?;

    // Call the cluster_routes function
    let clusters = cluster_routes(routes_data, threshold);

    // Serialize the clusters back to JsValue
    Ok(serde_wasm_bindgen::to_value(&clusters)?)
}

// Convert degrees to radians.
fn num2rad(d: &[f64]) -> Vec<f64> {
    d.iter().map(|&x| x.to_radians()).collect()
}

// Calculate the Haversine distance between two points in radians.
fn earth_haversine(p: &[f64], q: &[f64]) -> f64 {
    let d = [q[0] - p[0], q[1] - p[1]];
    let a = (d[0] / 2.0).sin().powi(2)
        + p[0].cos() * q[0].cos() * (d[1] / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
    c * 6371.0088 // Earth's radius in kilometers
}

// Calculate the Haversine distance between two points in degrees.
fn haversine_distance_deg(p: &[f64; 2], q: &[f64; 2]) -> f64 {
    let p_rad = num2rad(p);
    let q_rad = num2rad(q);
    earth_haversine(&p_rad, &q_rad)
}

// Compute the Frechet distance between two routes.
fn linear_frechet(p: &[[f64; 2]], q: &[[f64; 2]]) -> f64 {
    let num_p = p.len();
    let num_q = q.len();
    let mut ca = vec![vec![0.0; num_q]; num_p];

    for i in 0..num_p {
        for j in 0..num_q {
            let p_rad = num2rad(&p[i]);
            let q_rad = num2rad(&q[j]);
            let d = earth_haversine(&p_rad, &q_rad);

            if i > 0 && j > 0 {
                let min_prev = f64::min(
                    f64::min(ca[i - 1][j], ca[i - 1][j - 1]),
                    ca[i][j - 1],
                );
                ca[i][j] = f64::max(min_prev, d);
            } else if i > 0 && j == 0 {
                ca[i][j] = f64::max(ca[i - 1][0], d);
            } else if i == 0 && j > 0 {
                ca[i][j] = f64::max(ca[0][j - 1], d);
            } else {
                ca[i][j] = d;
            }
        }
    }
    ca[num_p - 1][num_q - 1]
}

// Compute the centroid of a route.
fn compute_centroid(points: &[[f64; 2]]) -> [f64; 2] {
    let n = points.len() as f64;
    let sum_lat = points.iter().map(|p| p[0]).sum::<f64>();
    let sum_lng = points.iter().map(|p| p[1]).sum::<f64>();
    [sum_lat / n, sum_lng / n]
}

// Get grid cell indices for spatial indexing.
fn get_grid_cell(
    lat: f64,
    lng: f64,
    min_lat: f64,
    min_lng: f64,
    cell_size: f64,
) -> (i32, i32) {
    let x = ((lat - min_lat) / cell_size).floor() as i32;
    let y = ((lng - min_lng) / cell_size).floor() as i32;
    (x, y)
}

// Cluster routes based on the Frechet distance.
fn cluster_routes(routes_data: Vec<Vec<[f64; 2]>>, threshold: f64) -> Vec<Vec<usize>> {
    let n = routes_data.len();
    let cell_size = 0.009; // Approximately 1 km in degrees
    let max_distance = 10.0; // Max distance between centroids in kilometers

    // Compute centroids and collect min_lat and min_lng
    let mut min_lat = f64::INFINITY;
    let mut min_lng = f64::INFINITY;

    let routes: Vec<_> = routes_data
        .into_iter()
        .enumerate()
        .map(|(i, points)| {
            let centroid = compute_centroid(&points);
            if centroid[0] < min_lat {
                min_lat = centroid[0];
            }
            if centroid[1] < min_lng {
                min_lng = centroid[1];
            }
            (i, points, centroid)
        })
        .collect();

    // Build grid map for spatial indexing
    let mut grid_map: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for (index, _points, centroid) in &routes {
        let (x, y) = get_grid_cell(centroid[0], centroid[1], min_lat, min_lng, cell_size);
        grid_map.entry((x, y)).or_default().push(*index);
    }

    let mut edges = Vec::new();

    // Process routes
    for (index, points, centroid) in &routes {
        let (x, y) = get_grid_cell(centroid[0], centroid[1], min_lat, min_lng, cell_size);

        // Check neighboring cells
        for dx in -1..=1 {
            for dy in -1..=1 {
                let neighbor_cell = (x + dx, y + dy);
                if let Some(neighbor_indices) = grid_map.get(&neighbor_cell) {
                    for &neighbor_index in neighbor_indices {
                        if neighbor_index > *index {
                            let neighbor = &routes[neighbor_index];
                            let d = haversine_distance_deg(centroid, &neighbor.2);
                            if d <= max_distance {
                                let frechet_d = linear_frechet(points, &neighbor.1);
                                if frechet_d <= threshold {
                                    edges.push((*index, neighbor_index));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Build the graph
    let mut graph = UnGraph::<usize, ()>::new_undirected();
    let node_indices: Vec<_> = (0..n).map(|i| graph.add_node(i)).collect();

    for &(i, j) in edges.iter() {
        graph.add_edge(node_indices[i], node_indices[j], ());
    }

    // Find connected components (clusters)
    let scc = kosaraju_scc(&graph);

    // Map NodeIndex to route indices
    scc.into_iter()
        .map(|component| component.into_iter().map(|node_index| graph[node_index]).collect())
        .collect()
}
