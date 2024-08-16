import * as sdf from "../SignedDistanceField.mjs";

// let model = sdf.cube({ position: [0,0,0], size: 20 });
// model = model.intersection_smooth(0.5, sdf.sphere({ position: [0, 0, 0], radius: 13 }));
// model = model.difference_smooth(0.5, sdf.cylinder({ position: [0, 0, 0], height: 20, radius: 6 }));
// model = model.difference_smooth(0.5, sdf.cylinder({ position: [0, 0, 0], height: 20, radius: 6 }).transform(sdf.Matrix4.rotatex(Math.PI / 2)));
// model = model.difference_smooth(0.5, sdf.cylinder({ position: [0, 0, 0], height: 20, radius: 6 }).transform(sdf.Matrix4.rotatey(Math.PI / 2)));
// console.log(model);
let model = sdf.cylinder({ position: [0, 0, 0], height: 20, radius: 6 }).transform(sdf.Matrix4.rotatex(Math.PI / 2));
console.log(model.calculate_signed_distance(50, 50, 50));
// const mesh = model.surface_nets(100);
// mesh.serialize_to_obj("example.obj");
