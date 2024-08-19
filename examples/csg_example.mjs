import { time, timeEnd } from "node:console";
import * as sdf from "../src/index.mjs";
import fs from "node:fs";

let model = sdf.cube({ center: [0,0,0], size: 20 });
model = model.intersection_smooth(0.5, sdf.sphere({ center: [0, 0, 0], radius: 13 }));
model = model.difference_smooth(0.5,
    sdf.cylinder({ center: [0, 0, 0], height: 20, radius: 6 }),
    sdf.cylinder({ center: [0, 0, 0], height: 20, radius: 6 }).rotate_x(Math.PI / 2),
    sdf.cylinder({ center: [0, 0, 0], height: 20, radius: 6 }).rotate_y(Math.PI / 2),
);

const mesh = model.surface_nets(200);
const buffer = mesh.serialize_to_stl();
console.log(buffer);
