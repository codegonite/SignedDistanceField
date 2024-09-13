import * as sdf from "../src/index.js";

const field = sdf.box(20, 20, 20)
    .intersectionSmooth(0.5, sdf.sphere(13))
    .differenceSmooth(0.5,
        sdf.cylinder(6, 20),
        sdf.cylinder(6, 20).rotateX(Math.PI / 2),
        sdf.cylinder(6, 20).rotateY(Math.PI / 2),
    );

const mesh = sdf.triangulateSignedDistanceField(field);
const buffer = sdf.convertToSTL(mesh);

console.log(buffer);
