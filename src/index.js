export const EPSILON = 1e-6;
export const INV_TWO_EPSILON = 1 / (2 * EPSILON);

const SAMPLE_OCTREE_KIND_FULL    = 0x01;
const SAMPLE_OCTREE_KIND_EMPTY   = 0x02;
const SAMPLE_OCTREE_KIND_ROOT    = 0x04;
const SAMPLE_OCTREE_KIND_LEAF    = 0x08;

const SAMPLE_OCTREE_FULL_OR_EMPTY = SAMPLE_OCTREE_KIND_FULL | SAMPLE_OCTREE_KIND_EMPTY;

const ALL_SAMPLES_INSIDE = 0xFF;
const ALL_SAMPLES_OUTSIDE = 0x00;

const STL_HEADER_SIZE = 80;
const STL_TRIANGLE_OFFSET = 84;
const STL_BYTES_PER_TRIANGLE = 50;

const INTERSECTION_DEPTH_MAX = 256;

export class Vector2 {
    static fromObject(v) { return new Vector2(v.x, v.y); }
    static fromArray(v)  { return new Vector2(v[0], v[1]); }

    constructor(x = 0, y = 0) {
        this.x = x;
        this.y = y;
    }

    *[Symbol.iterator]() {
        yield this.x;
        yield this.y;
    }

    clone() { return new Vector2(this.x, this.y); }
    toObject() { return { x: this.x, y: this.y }; }
    toArray() { return [this.x, this.y]; }

    distance(other) {
        const dx = other.x - this.x;
        const dy = other.y - this.y;
        return Math.sqrt(dx * dx + dy * dy);
    }

    magnitudeSquared() {
        return this.x * this.x + this.y * this.y;
    }

    magnitude() {
        return Math.sqrt(this.x * this.x + this.y * this.y);
    }

    dot(other) {
        return this.x * other.x + this.y * other.y;
    }

    apply(callback) {
        this.x = callback(this.x);
        this.y = callback(this.y);
        return this;
    }

    min(other) {
        this.x = Math.min(this.x, other.x);
        this.y = Math.min(this.y, other.y);
        return this;
    }

    max(other) {
        this.x = Math.max(this.x, other.x);
        this.y = Math.max(this.y, other.y);
        return this;
    }

    clamp(v0, v1) {
        this.x = clamp(this.x, v0.x, v1.x);
        this.y = clamp(this.y, v0.y, v1.y);
        return this;
    }

    neg() {
        this.x = -this.x;
        this.y = -this.y;
        return this;
    }

    add(other) {
        this.x += other.x;
        this.y += other.y;
        return this;
    }

    sub(other) {
        this.x -= other.x;
        this.y -= other.y;
        return this;
    }

    mul(other) {
        this.x *= other.x;
        this.y *= other.y;
        return this;
    }

    div(other) {
        this.x /= other.x;
        this.y /= other.y;
        return this;
    }

    addScalar(scalar) {
        this.x += scalar;
        this.y += scalar;
        return this;
    }

    subScalar(scalar) {
        this.x -= scalar;
        this.y -= scalar;
        return this;
    }

    mulScalar(scalar) {
        this.x *= scalar;
        this.y *= scalar;
        return this;
    }

    divScalar(scalar) {
        this.x /= scalar;
        this.y /= scalar;
        return this;
    }

    setLength(length) {
        const length_squared = this.x * this.x + this.y * this.y;
        if (length_squared === 0) {
            this.x = this.y = 0;
            return this;
        }
        const factor = length / Math.sqrt(length_squared);
        this.x *= factor;
        this.y *= factor;
        return this;
    }

    normalize() {
        const length_squared = this.x * this.x + this.y * this.y;
        if (length_squared === 0) {
            this.x = this.y = 0;
            return this;
        }
        const factor = 1 / Math.sqrt(length_squared);
        this.x *= factor;
        this.y *= factor;
        return this;
    }

    transformMatrix4(matrix) {
        const x = this.x, y = this.y;
        this.x = x * matrix.m00 + y * matrix.m10 + matrix.m20;
        this.y = x * matrix.m01 + y * matrix.m11 + matrix.m21;
        return this;
    }

    transformMatrix3(matrix) {
        const x = this.x, y = this.y;
        this.x = x * matrix.m00 + y * matrix.m10 + matrix.m20;
        this.y = x * matrix.m01 + y * matrix.m11 + matrix.m21;
        return this;
    }

    transformMatrix2(matrix) {
        const x = this.x, y = this.y;
        this.x = x * matrix.m00 + y * matrix.m10;
        this.y = x * matrix.m01 + y * matrix.m11;
        return this;
    }

    cross() {
        return this.x * this.y - this.y * this.x;
    }

    set(x, y) {
        this.x = x;
        this.y = y;
        return this;
    }

    copy(other) {
        this.x = other.x;
        this.y = other.y;
        return this;
    }

    addScaledVector(other, scalar) {
        this.x += other.x * scalar;
        this.y += other.y * scalar;
        return this;
    }

    subScaledVector(other, scalar) {
        this.x -= other.x * scalar;
        this.y -= other.y * scalar;
        return this;
    }

    negVector(other) {
        this.x = -other.x;
        this.y = -other.y;
        return this;
    }

    addVectors(v0, v1) {
        this.x = v0.x + v1.x;
        this.y = v0.y + v1.y;
        return this;
    }

    subVectors(v0, v1) {
        this.x = v0.x - v1.x;
        this.y = v0.y - v1.y;
        return this;
    }

    mulVectors(v0, v1) {
        this.x = v0.x * v1.x;
        this.y = v0.y * v1.y;
        return this;
    }

    divVectors(v0, v1) {
        this.x = v0.x / v1.x;
        this.y = v0.y / v1.y;
        return this;
    }

    lerpVectors(v0, v1, t = 0) {
        this.x = (v1.x - v0.x) * t + v0.x;
        this.y = (v1.y - v0.y) * t + v0.y;
        return this;
    }
}

export class Vector3 {
    static fromObject(v) { return new Vector3(v.x, v.y, v.z); }
    static fromArray(v)  { return new Vector3(v[0], v[1], v[2]); }

    constructor(x=0, y=0, z=0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    *[Symbol.iterator]() {
        yield this.x;
        yield this.y;
        yield this.z;
    }

    clone() { return new Vector3(this.x, this.y, this.z); }
    toObject() { return { x: this.x, y: this.y, z: this.z }; }
    toArray() { return [ this.x, this.y, this.z ]; }

    distance(other) {
        const x = other.x - this.x;
        const y = other.y - this.y;
        const z = other.z - this.z;
        return Math.sqrt(x*x + y*y + z*z);
    }

    magnitudeSquared() {
        return this.x*this.x + this.y*this.y + this.z*this.z;
    }
    
    magnitude() {
        return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    }
    
    dot(other) {
        return this.x*other.x + this.y*other.y + this.z*other.z;
    }

    apply(callback) {
        this.x = callback(this.x);
        this.y = callback(this.y);
        this.z = callback(this.z);
        return this;
    }

    min(other) {
        this.x = Math.min(this.x, other.x);
        this.y = Math.min(this.y, other.y);
        this.z = Math.min(this.z, other.z);
        return this;
    }

    max(other) {
        this.x = Math.max(this.x, other.x);
        this.y = Math.max(this.y, other.y);
        this.z = Math.max(this.z, other.z);
        return this;
    }

    clamp(v0, v1) {
        this.x = clamp(this.x, v0.x, v1.x);
        this.y = clamp(this.y, v0.y, v1.y);
        this.z = clamp(this.z, v0.z, v1.z);
        return this;
    }

    neg() {
        this.x = -this.x;
        this.y = -this.y;
        this.z = -this.z;
        return this;
    }

    add(other) {
        this.x += other.x;
        this.y += other.y;
        this.z += other.z;
        return this;
    }

    sub(other) {
        this.x -= other.x;
        this.y -= other.y;
        this.z -= other.z;
        return this;
    }

    mul(other) {
        this.x *= other.x;
        this.y *= other.y;
        this.z *= other.z;
        return this;
    }

    div(other) {
        this.x /= other.x;
        this.y /= other.y;
        this.z /= other.z;
        return this;
    }

    addScalar(scalar) {
        this.x += scalar;
        this.y += scalar;
        this.z += scalar;
        return this;
    }

    subScalar(scalar) {
        this.x -= scalar;
        this.y -= scalar;
        this.z -= scalar;
        return this;
    }

    mulScalar(scalar) {
        this.x *= scalar;
        this.y *= scalar;
        this.z *= scalar;
        return this;
    }

    divScalar(scalar) {
        this.x /= scalar;
        this.y /= scalar;
        this.z /= scalar;
        return this;
    }

    setLength(length) {
        const x0=this.x, y0=this.y, z0=this.z;
        const length_squared = x0*x0 + y0*y0 + z0*z0;

        if (length_squared == 0) {
            this.x = this.y = this.z = 0;
            return this;
        }

        const factor = length / Math.sqrt(length_squared);
        this.x = x0 * factor;
        this.y = y0 * factor;
        this.z = z0 * factor;
        return this;
    }

    normalize() {
        const x0=this.x, y0=this.y, z0=this.z;
        const length_squared = x0*x0 + y0*y0 + z0*z0;

        if (length_squared !== 0) {
            const factor = 1 / Math.sqrt(length_squared);
            this.x = x0 * factor;
            this.y = y0 * factor;
            this.z = z0 * factor;
            return this;
        }

        this.x = this.y = this.z = 0;
        return this;
    }

    transformMatrix3(m) {
        const x = this.x, y = this.y, z = this.z;
        this.x = x*m.m00 + y*m.m10 + z*m.m20;
        this.y = x*m.m01 + y*m.m11 + z*m.m21;
        this.z = x*m.m02 + y*m.m12 + z*m.m22;
        return this;
    }

    transformMatrix4(m) {
        const x = this.x, y = this.y, z = this.z;
        this.x = x*m.m00 + y*m.m10 + z*m.m20 + m.m30;
        this.y = x*m.m01 + y*m.m11 + z*m.m21 + m.m31;
        this.z = x*m.m02 + y*m.m12 + z*m.m22 + m.m32;
        return this;
    }

    cross(other) {
        const x0 = this.x,  y0 = this.y,  z0 = this.z;
        const x1 = other.x, y1 = other.y, z1 = other.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }

    set(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }

    copy(other) {
        this.x = other.x;
        this.y = other.y;
        this.z = other.z;
        return this;
    }

    addScaledVector(other, scalar) {
        this.x += other.x * scalar;
        this.y += other.y * scalar;
        this.z += other.z * scalar;
        return this;
    }

    subScaledVector(other, scalar) {
        this.x -= other.x * scalar;
        this.y -= other.y * scalar;
        this.z -= other.z * scalar;
        return this;
    }


    negVector(other) {
        this.x = -other.x;
        this.y = -other.y;
        this.z = -other.z;
        return this;
    }

    addVectors(v0, v1) {
        this.x = v0.x+v1.x;
        this.y = v0.y+v1.y;
        this.z = v0.z+v1.z;
        return this;
    }

    subVectors(v0, v1) {
        this.x = v0.x-v1.x;
        this.y = v0.y-v1.y;
        this.z = v0.z-v1.z;
        return this;
    }

    mulVectors(v0, v1) {
        this.x = v0.x*v1.x;
        this.y = v0.y*v1.y;
        this.z = v0.z*v1.z;
        return this;
    }

    divVectors(v0, v1) {
        this.x = v0.x/v1.x;
        this.y = v0.y/v1.y;
        this.z = v0.z/v1.z;
        return this;
    }

    lerpVectors(v0, v1, t = 0) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;

        this.x = (x1 - x0) * t + x0;
        this.y = (y1 - y0) * t + y0;
        this.z = (z1 - z0) * t + z0;
        return this;
    }

    crossVectors(v0, v1) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }
}

export class Vector4 {
    static fromObject(v) { return new Vector4(v.x, v.y, v.z, v.w); }
    static fromArray(v)  { return new Vector4(v[0], v[1], v[2], v[3]); }

    constructor(x = 0, y = 0, z = 0, w = 1) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    *[Symbol.iterator]() {
        yield this.x;
        yield this.y;
        yield this.z;
        yield this.w;
    }

    clone() { return new Vector4(this.x, this.y, this.z, this.w); }
    toObject() { return { x: this.x, y: this.y, z: this.z, w: this.w }; }
    toArray() { return [this.x, this.y, this.z, this.w]; }

    distance(other) {
        const dx = other.x - this.x;
        const dy = other.y - this.y;
        const dz = other.z - this.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    magnitudeSquared() {
        return this.x * this.x + this.y * this.y + this.z * this.z;
    }

    magnitude() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }

    dot(other) {
        return this.x * other.x + this.y * other.y + this.z * other.z + this.w * other.w;
    }

    apply(callback) {
        this.x = callback(this.x);
        this.y = callback(this.y);
        this.z = callback(this.z);
        this.w = callback(this.w);
        return this;
    }

    min(other) {
        this.x = Math.min(this.x, other.x);
        this.y = Math.min(this.y, other.y);
        this.z = Math.min(this.z, other.z);
        this.w = Math.min(this.w, other.w);
        return this;
    }

    max(other) {
        this.x = Math.max(this.x, other.x);
        this.y = Math.max(this.y, other.y);
        this.z = Math.max(this.z, other.z);
        this.w = Math.max(this.w, other.w);
        return this;
    }

    clamp(v0, v1) {
        this.x = clamp(this.x, v0.x, v1.x);
        this.y = clamp(this.y, v0.y, v1.y);
        this.z = clamp(this.z, v0.z, v1.z);
        this.w = clamp(this.w, v0.w, v1.w);
        return this;
    }

    neg() {
        this.x = -this.x;
        this.y = -this.y;
        this.z = -this.z;
        this.w = -this.w;
        return this;
    }

    add(other) {
        this.x += other.x;
        this.y += other.y;
        this.z += other.z;
        this.w += other.w;
        return this;
    }

    sub(other) {
        this.x -= other.x;
        this.y -= other.y;
        this.z -= other.z;
        this.w -= other.w;
        return this;
    }

    mul(other) {
        this.x *= other.x;
        this.y *= other.y;
        this.z *= other.z;
        this.w *= other.w;
        return this;
    }

    div(other) {
        this.x /= other.x;
        this.y /= other.y;
        this.z /= other.z;
        this.w /= other.w;
        return this;
    }

    addScalar(scalar) {
        this.x += scalar;
        this.y += scalar;
        this.z += scalar;
        this.w += scalar;
        return this;
    }

    subScalar(scalar) {
        this.x -= scalar;
        this.y -= scalar;
        this.z -= scalar;
        this.w -= scalar;
        return this;
    }

    mulScalar(scalar) {
        this.x *= scalar;
        this.y *= scalar;
        this.z *= scalar;
        this.w *= scalar;
        return this;
    }

    divScalar(scalar) {
        this.x /= scalar;
        this.y /= scalar;
        this.z /= scalar;
        this.w /= scalar;
        return this;
    }

    setLength(length) {
        const length_squared = this.x * this.x + this.y * this.y + this.z * this.z;
        if (length_squared === 0) {
            this.x = 0;
            this.y = 0;
            this.z = 0;
            this.w = 1;
            return this;
        }
        const factor = length / Math.sqrt(length_squared);
        this.x *= factor;
        this.y *= factor;
        this.z *= factor;
        this.w *= factor;
        return this;
    }

    normalize() {
        const length_squared = this.x * this.x + this.y * this.y + this.z * this.z;
        if (length_squared === 0) {
            this.x = 0;
            this.y = 0;
            this.z = 0;
            this.w = 1;
            return this;
        }
        const factor = 1 / Math.sqrt(length_squared);
        this.x *= factor;
        this.y *= factor;
        this.z *= factor;
        this.w *= factor;
        return this;
    }

    transformMatrix4(m) {
        const x = this.x, y = this.y, z = this.z, w = this.w;
        this.x = x * m.m00 + y * m.m10 + z * m.m20 + w * m.m30;
        this.y = x * m.m01 + y * m.m11 + z * m.m21 + w * m.m31;
        this.z = x * m.m02 + y * m.m12 + z * m.m22 + w * m.m32;
        this.w = x * m.m03 + y * m.m13 + z * m.m23 + w * m.m33;
        return this;
    }

    cross(other) {
        const x0 = this.x, y0 = this.y, z0 = this.z, w0 = this.w;
        const x1 = other.x, y1 = other.y, z1 = other.z, w1 = other.w;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        this.w = w0;
        return this;
    }

    set(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
        return this;
    }

    copy(other) {
        this.x = other.x;
        this.y = other.y;
        this.z = other.z;
        this.w = other.w;
        return this;
    }

    addScaledVector(other, scalar) {
        this.x += other.x * scalar;
        this.y += other.y * scalar;
        this.z += other.z * scalar;
        this.w += other.w * scalar;
        return this;
    }

    subScaledVector(other, scalar) {
        this.x -= other.x * scalar;
        this.y -= other.y * scalar;
        this.z -= other.z * scalar;
        this.w -= other.w * scalar;
        return this;
    }

    negVector(other) {
        this.x = -other.x;
        this.y = -other.y;
        this.z = -other.z;
        this.w = -other.w;
        return this;
    }

    addVectors(v0, v1) {
        this.x = v0.x + v1.x;
        this.y = v0.y + v1.y;
        this.z = v0.z + v1.z;
        this.w = v0.w + v1.w;
        return this;
    }

    subVectors(v0, v1) {
        this.x = v0.x - v1.x;
        this.y = v0.y - v1.y;
        this.z = v0.z - v1.z;
        this.w = v0.w - v1.w;
        return this;
    }

    mulVectors(v0, v1) {
        this.x = v0.x * v1.x;
        this.y = v0.y * v1.y;
        this.z = v0.z * v1.z;
        this.w = v0.w * v1.w;
        return this;
    }

    divVectors(v0, v1) {
        this.x = v0.x / v1.x;
        this.y = v0.y / v1.y;
        this.z = v0.z / v1.z;
        this.w = v0.w / v1.w;
        return this;
    }

    lerpVectors(v0, v1, t = 0) {
        this.x = (v1.x - v0.x) * t + v0.x;
        this.y = (v1.y - v0.y) * t + v0.y;
        this.z = (v1.z - v0.z) * t + v0.z;
        this.w = (v1.w - v0.w) * t + v0.w;
        return this;
    }

    crossVectors(v0, v1) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z, w0 = v0.w;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z, w1 = v1.w;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        this.w = w0;
        return this;
    }
}

export class BoundingBox2 {
    constructor(min_x = 0, min_y = 0, max_x = 0, max_y = 0) {
        this.min_x = min_x;
        this.max_x = max_x;
        this.min_y = min_y;
        this.max_y = max_y;
    }

    clone() {
        return new BoundingBox2(this.min_x, this.min_y, this.max_x, this.max_y);
    }

    set(min_x = 0, min_y = 0, max_x = 0, max_y = 0) {
        this.min_x = min_x;
        this.max_x = max_x;
        this.min_y = min_y;
        this.max_y = max_y;
        return this;
    }

    copy(box) {
        this.min_x = box.min_x;
        this.max_x = box.max_x;
        this.min_y = box.min_y;
        this.max_y = box.max_y;
        return this;
    }

    grow(size) {
        this.min_x -= size;
        this.min_y -= size;
        this.max_x += size;
        this.max_y += size;
        return this;
    }

    mul(factor) {
        this.min_x *= factor;
        this.min_y *= factor;
        this.max_x *= factor;
        this.max_y *= factor;
        return this;
    }

    transformMatrix2(matrix) {
        const r0 = matrix.m00, r1 = matrix.m01;
        const u0 = matrix.m10, u1 = matrix.m11;

        const xa0 = r0 * this.min_x, xa1 = r1 * this.min_x;
        const xb0 = r0 * this.max_x, xb1 = r1 * this.max_x;
        const ya0 = u0 * this.min_y, ya1 = u1 * this.min_y;
        const yb0 = u0 * this.max_y, yb1 = u1 * this.max_y;

        const min0_x = Math.min(xa0, xb0), min1_x = Math.min(xa1, xb1);
        const max0_x = Math.max(xa0, xb0), max1_x = Math.max(xa1, xb1);
        const min0_y = Math.min(ya0, yb0), min1_y = Math.min(ya1, yb1);
        const max0_y = Math.max(ya0, yb0), max1_y = Math.max(ya1, yb1);

        this.min_x = min0_x + min0_y;
        this.min_y = min1_x + min1_y;
        this.max_x = max0_x + max0_y;
        this.max_y = max1_x + max1_y;
        return this;
    }

    transformMatrix3(matrix) {
        const r0 = matrix.m00, r1 = matrix.m01;
        const u0 = matrix.m10, u1 = matrix.m11;
        const t0 = matrix.m20, t1 = matrix.m21;

        const xa0 = r0 * this.min_x, xa1 = r1 * this.min_x;
        const xb0 = r0 * this.max_x, xb1 = r1 * this.max_x;
        const ya0 = u0 * this.min_y, ya1 = u1 * this.min_y;
        const yb0 = u0 * this.max_y, yb1 = u1 * this.max_y;

        const min0_x = Math.min(xa0, xb0), min1_x = Math.min(xa1, xb1);
        const max0_x = Math.max(xa0, xb0), max1_x = Math.max(xa1, xb1);
        const min0_y = Math.min(ya0, yb0), min1_y = Math.min(ya1, yb1);
        const max0_y = Math.max(ya0, yb0), max1_y = Math.max(ya1, yb1);

        this.min_x = min0_x + min0_y + t0;
        this.min_y = min1_x + min1_y + t1;
        this.max_x = max0_x + max0_y + t0;
        this.max_y = max1_x + max1_y + t1;
        return this;
    }

    unionBoxes(box0, box1) {
        this.min_x = Math.min(box0.min_x, box1.min_x);
        this.min_y = Math.min(box0.min_y, box1.min_y);
        this.max_x = Math.max(box0.max_x, box1.max_x);
        this.max_y = Math.max(box0.max_y, box1.max_y);
        return this;
    }

    intersectBoxes(box0, box1) {
        this.min_x = Math.max(box0.min_x, box1.min_x);
        this.min_y = Math.max(box0.min_y, box1.min_y);
        this.max_x = Math.min(box0.max_x, box1.max_x);
        this.max_y = Math.min(box0.max_y, box1.max_y);
        return this;
    }

    union(box) {
        return this.unionBoxes(this, box);
    }

    intersect(box) {
        return this.intersectBoxes(this, box);
    }
}

export class BoundingBox3 {
    constructor(min_x = 0, min_y = 0, min_z = 0, max_x = 0, max_y = 0, max_z = 0) {
        this.min_x = min_x;
        this.min_y = min_y;
        this.min_z = min_z;
        this.max_x = max_x;
        this.max_y = max_y;
        this.max_z = max_z;
    }

    clone() {
        return new BoundingBox3(this.min_x, this.min_y, this.min_z,
                                this.max_x, this.max_y, this.max_z);
    }

    set(min_x = 0, min_y = 0, min_z = 0, max_x = 0, max_y = 0, max_z = 0) {
        this.min_x = min_x;
        this.min_y = min_y;
        this.min_z = min_z;
        this.max_x = max_x;
        this.max_y = max_y;
        this.max_z = max_z;
        return this;
    }

    copy(box) {
        this.min_x = box.min_x;
        this.min_y = box.min_y;
        this.min_z = box.min_z;
        this.max_x = box.max_x;
        this.max_y = box.max_y;
        this.max_z = box.max_z;
        return this;
    }

    grow(size) {
        this.min_x -= size;
        this.min_y -= size;
        this.min_z -= size;
        this.max_x += size;
        this.max_y += size;
        this.max_z += size;
        return this;
    }

    mul(factor) {
        this.min_x *= factor;
        this.min_y *= factor;
        this.min_z *= factor;
        this.max_x *= factor;
        this.max_y *= factor;
        this.max_z *= factor;
        return this;
    }

    transformMatrix3(matrix) {
        const r0 = matrix.m00, r1 = matrix.m01, r2 = matrix.m02;
        const u0 = matrix.m10, u1 = matrix.m11, u2 = matrix.m12;
        const b0 = matrix.m20, b1 = matrix.m21, b2 = matrix.m22;

        const xa0 = r0 * this.min_x, xa1 = r1 * this.min_x, xa2 = r2 * this.min_x;
        const xb0 = r0 * this.max_x, xb1 = r1 * this.max_x, xb2 = r2 * this.max_x;
        const ya0 = u0 * this.min_y, ya1 = u1 * this.min_y, ya2 = u2 * this.min_y;
        const yb0 = u0 * this.max_y, yb1 = u1 * this.max_y, yb2 = u2 * this.max_y;
        const za0 = b0 * this.min_z, za1 = b1 * this.min_z, za2 = b2 * this.min_z;
        const zb0 = b0 * this.max_z, zb1 = b1 * this.max_z, zb2 = b2 * this.max_z;

        const min0_x = Math.min(xa0, xb0), min1_x = Math.min(xa1, xb1), min2_x = Math.min(xa2, xb2);
        const max0_x = Math.max(xa0, xb0), max1_x = Math.max(xa1, xb1), max2_x = Math.max(xa2, xb2);
        const min0_y = Math.min(ya0, yb0), min1_y = Math.min(ya1, yb1), min2_y = Math.min(ya2, yb2);
        const max0_y = Math.max(ya0, yb0), max1_y = Math.max(ya1, yb1), max2_y = Math.max(ya2, yb2);
        const min0_z = Math.min(za0, zb0), min1_z = Math.min(za1, zb1), min2_z = Math.min(za2, zb2);
        const max0_z = Math.max(za0, zb0), max1_z = Math.max(za1, zb1), max2_z = Math.max(za2, zb2);

        this.min_x = min0_x + min0_y + min0_z;
        this.min_y = min1_x + min1_y + min1_z;
        this.min_z = min2_x + min2_y + min2_z;
        this.max_x = max0_x + max0_y + max0_z;
        this.max_y = max1_x + max1_y + max1_z;
        this.max_z = max2_x + max2_y + max2_z;
        return this;
    }
    
    transformMatrix4(matrix) {
        const r0 = matrix.m00, r1 = matrix.m01, r2 = matrix.m02;
        const u0 = matrix.m10, u1 = matrix.m11, u2 = matrix.m12;
        const b0 = matrix.m20, b1 = matrix.m21, b2 = matrix.m22;
        const t0 = matrix.m30, t1 = matrix.m31, t2 = matrix.m32;

        const xa0 = r0 * this.min_x, xa1 = r1 * this.min_x, xa2 = r2 * this.min_x;
        const xb0 = r0 * this.max_x, xb1 = r1 * this.max_x, xb2 = r2 * this.max_x;
        const ya0 = u0 * this.min_y, ya1 = u1 * this.min_y, ya2 = u2 * this.min_y;
        const yb0 = u0 * this.max_y, yb1 = u1 * this.max_y, yb2 = u2 * this.max_y;
        const za0 = b0 * this.min_z, za1 = b1 * this.min_z, za2 = b2 * this.min_z;
        const zb0 = b0 * this.max_z, zb1 = b1 * this.max_z, zb2 = b2 * this.max_z;

        const min0_x = Math.min(xa0, xb0), min1_x = Math.min(xa1, xb1), min2_x = Math.min(xa2, xb2);
        const max0_x = Math.max(xa0, xb0), max1_x = Math.max(xa1, xb1), max2_x = Math.max(xa2, xb2);
        const min0_y = Math.min(ya0, yb0), min1_y = Math.min(ya1, yb1), min2_y = Math.min(ya2, yb2);
        const max0_y = Math.max(ya0, yb0), max1_y = Math.max(ya1, yb1), max2_y = Math.max(ya2, yb2);
        const min0_z = Math.min(za0, zb0), min1_z = Math.min(za1, zb1), min2_z = Math.min(za2, zb2);
        const max0_z = Math.max(za0, zb0), max1_z = Math.max(za1, zb1), max2_z = Math.max(za2, zb2);

        this.min_x = min0_x + min0_y + min0_z + t0;
        this.min_y = min1_x + min1_y + min1_z + t1;
        this.min_z = min2_x + min2_y + min2_z + t2;
        this.max_x = max0_x + max0_y + max0_z + t0;
        this.max_y = max1_x + max1_y + max1_z + t1;
        this.max_z = max2_x + max2_y + max2_z + t2;
        return this;
    }

    unionBoxes(box0, box1) {
        this.min_x = Math.min(box0.min_x, box1.min_x);
        this.min_y = Math.min(box0.min_y, box1.min_y);
        this.min_z = Math.min(box0.min_z, box1.min_z);
        this.max_x = Math.max(box0.max_x, box1.max_x);
        this.max_y = Math.max(box0.max_y, box1.max_y);
        this.max_z = Math.max(box0.max_z, box1.max_z);
        return this;
    }

    intersectBoxes(box0, box1) {
        this.min_x = Math.max(box0.min_x, box1.min_x);
        this.min_y = Math.max(box0.min_y, box1.min_y);
        this.min_z = Math.max(box0.min_z, box1.min_z);
        this.max_x = Math.min(box0.max_x, box1.max_x);
        this.max_y = Math.min(box0.max_y, box1.max_y);
        this.max_z = Math.min(box0.max_z, box1.max_z);
        return this;
    }

    union(box) {
        return this.unionBoxes(this, box);
    }

    intersect(box) {
        return this.intersectBoxes(this, box);
    }
}

export class Plane {
    constructor(normal_x, normal_y, normal_z, distance) {
        this.normal_x = normal_x;
        this.normal_y = normal_y;
        this.normal_z = normal_z;
        this.distance = distance;
    }

    getOrigin(out = new Vector3()) {
        return out.set(this.normal_x * this.distance,
                       this.normal_y * this.distance,
                       this.normal_z * this.distance);
    }

    getNormal(out = new Vector3()) {
        return out.set(this.normal_x,
                       this.normal_y,
                       this.normal_z);
    }

    clone() {
        return new Plane(this.normal_x, this.normal_y, this.normal_z, this.distance);
    }

    dotPlane(plane) {
        return this.normal_x * plane.normal_x + this.normal_y * plane.normal_y + this.normal_z * plane.normal_z;
    }

    dotVector(vector) {
        return this.normal_x * vector.x + this.normal_y * vector.y + this.normal_z * vector.z;
    }

    copy(other) {
        this.normal_x = other.normal_x;
        this.normal_y = other.normal_y;
        this.normal_z = other.normal_z;
        this.distance = other.distance;
        return this;
    }

    equals(other) {
        const x = other.normal_x - this.normal_x;
        const y = other.normal_y - this.normal_y;
        const z = other.normal_z - this.normal_z;
        const plane_distance = Math.abs(this.distance - other.distance);
        const origin_distance = Math.sqrt(x*x + y*y + z*z);
        return plane_distance <= EPSILON && origin_distance <= EPSILON;
    }

    flip() {
        this.normal_x = -this.normal_x;
        this.normal_y = -this.normal_y;
        this.normal_z = -this.normal_z;
        this.distance = -this.distance;
        return this;
    }

    distanceToPoint(x, y, z) {
        return x * this.normal_x + y * this.normal_y + z * this.normal_z - this.distance;
    }

    projectionOfPoint(x, y, z) {
        const a = x * this.normal_x + y * this.normal_y + z * this.normal_z - this.distance;
        return new Vector3(x - a * this.normal_x,
                           y - a * this.normal_y,
                           z - a * this.normal_z);
    }

    transformMatrix3(matrix, inverse_matrix) {
        const x=this.normal_x, y=this.normal_y, z=this.normal_z, d=this.distance;

        const origin0 = x*d*matrix.m00 + y*d*matrix.m10 + z*d*matrix.m20;
        const origin1 = x*d*matrix.m01 + y*d*matrix.m11 + z*d*matrix.m21;
        const origin2 = x*d*matrix.m02 + y*d*matrix.m12 + z*d*matrix.m22;

        const normal0 = x*inverse_matrix.m00 + y*inverse_matrix.m10 + z*inverse_matrix.m20;
        const normal1 = x*inverse_matrix.m01 + y*inverse_matrix.m11 + z*inverse_matrix.m21;
        const normal2 = x*inverse_matrix.m02 + y*inverse_matrix.m12 + z*inverse_matrix.m22;

        const distance = normal0*origin0 + normal1*origin1 + normal2*origin2;

        this.normal_x = normal0;
        this.normal_y = normal1;
        this.normal_z = normal2;
        this.distance = distance;
        return this;
    }

    transformMatrix4(matrix, inverse_matrix) {
        const x=this.normal_x, y=this.normal_y, z=this.normal_z, d=this.distance;

        const origin0 = x*d*matrix.m00 + y*d*matrix.m10 + z*d*matrix.m20 + matrix.m30;
        const origin1 = x*d*matrix.m01 + y*d*matrix.m11 + z*d*matrix.m21 + matrix.m31;
        const origin2 = x*d*matrix.m02 + y*d*matrix.m12 + z*d*matrix.m22 + matrix.m32;

        const normal0 = x*inverse_matrix.m00 + y*inverse_matrix.m10 + z*inverse_matrix.m20;
        const normal1 = x*inverse_matrix.m01 + y*inverse_matrix.m11 + z*inverse_matrix.m21;
        const normal2 = x*inverse_matrix.m02 + y*inverse_matrix.m12 + z*inverse_matrix.m22;

        const distance = normal0*origin0 + normal1*origin1 + normal2*origin2;

        this.normal_x = normal0;
        this.normal_y = normal1;
        this.normal_z = normal2;
        this.distance = distance;
        return this;
    }

    fromNormalAndPoint(normal_x, normal_y, normal_z, point_x, point_y, point_z) {
        const distance = normal_x*point_x + normal_y*point_y + normal_z*point_z;
        this.normal_x = normal_x;
        this.normal_y = normal_y;
        this.normal_z = normal_z;
        this.distance = distance;
        return this;
    }

    fromPoints(x0, y0, z0, x1, y1, z1, x2, y2, z2) {
        const dx0     = x1 - x0;
        const dy0     = y1 - y0;
        const dz0     = z1 - z0;
        const dx1     = x2 - x0;
        const dy1     = y2 - y0;
        const dz1     = z2 - z0;
        const cross_x = dy0 * dz1 - dy1 * dz0;
        const cross_y = dz0 * dx1 - dz1 * dx0;
        const cross_z = dx0 * dy1 - dx1 * dy0;
        const length_squared = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;

        if (length_squared === 0) {
            return new Plane(0, 0, 0, 0);
        }

        const factor   = 1.0 / Math.sqrt(length_squared);
        const normal_x = cross_x * factor;
        const normal_y = cross_y * factor;
        const normal_z = cross_z * factor;
        const distance = normal_x*x0 + normal_y*y0 + normal_z*z0;

        this.normal_x = normal_x;
        this.normal_y = normal_y;
        this.normal_z = normal_z;
        this.distance = distance;
        return this;
    }

    fromVectors(v0, v1, v2) {
        return this.fromPoints(v0.x, v0.y, v0.z,
                               v1.x, v1.y, v1.z,
                               v2.x, v2.y, v2.z);
    }
}

export class Matrix2 {
    static identity() {
        return new Matrix2(1, 0,
                           0, 1);
    }

    static scale(x, y) {
        return new Matrix2(x, 0,
                           0, y);
    }

    static rotate(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix2(c, -s,
                           s,  c);
    }

    constructor(m00=1, m10=0, m01=0, m11=1) {
        this.m00 = m00;
        this.m01 = m01;
        this.m10 = m10;
        this.m11 = m11;
    }

    clone() {
        return new Matrix2(this.m00, this.m10,
                           this.m01, this.m11);
    }

    determinant() {
        return this.m00 * this.m11 - this.m10 * this.m01;
    }

    isIdentity() {
        return this.m00 === 1 && this.m01 === 0
        &&     this.m10 === 0 && this.m11 === 1;
    }

    equals(other) {
        return Math.abs(this.m00 - other.m00) <= EPSILON && Math.abs(this.m10 - other.m10) <= EPSILON
        &&     Math.abs(this.m01 - other.m01) <= EPSILON && Math.abs(this.m11 - other.m11) <= EPSILON;
    }
    
    identity() {
        this.m00 = 1;
        this.m01 = 0;
        this.m10 = 0;
        this.m11 = 1;
        return this;
    }

    set(m00, m01, m10, m11) {
        this.m00 = m00;
        this.m01 = m01;
        this.m10 = m10;
        this.m11 = m11;
        return this;
    }

    copy(other) {
        this.m00 = other.m00;
        this.m01 = other.m01;
        this.m10 = other.m10;
        this.m11 = other.m11;
        return this;
    }

    transposeMatrix(matrix) {
        const m00 = matrix.m00, m01 = matrix.m01;
        const m10 = matrix.m10, m11 = matrix.m11;

        this.m00 = m00;
        this.m01 = m10;
        this.m10 = m01;
        this.m11 = m11;
        return this;
    }

    mulMatrices(matrix0, matrix1) {
        const a00 = matrix0.m00, a01 = matrix0.m01;
        const a10 = matrix0.m10, a11 = matrix0.m11;

        const b00 = matrix1.m00, b01 = matrix1.m01;
        const b10 = matrix1.m10, b11 = matrix1.m11;

        this.m00 = a00 * b00 + a01 * b10;
        this.m01 = a00 * b01 + a01 * b11;
        this.m10 = a10 * b00 + a11 * b10;
        this.m11 = a10 * b01 + a11 * b11;
        return this;
    }

    inverseMatrix(matrix) {
        const m00 = matrix.m00, m01 = matrix.m01;
        const m10 = matrix.m10, m11 = matrix.m11;

        const determinant = m00 * m11 - m10 * m01;

        if (determinant === 0) {
            throw new Error("Determinant of the matrix cannot be zero!");
        }

        const factor = 1 / determinant;

        this.m00 = m11 * factor;
        this.m01 = m10 * factor;
        this.m10 = m01 * factor;
        this.m11 = m00 * factor;

        return this;
    }

    scaleMatrix(matrix, x, y) {
        const m00 = matrix.m00, m01 = matrix.m01;
        const m10 = matrix.m10, m11 = matrix.m11;

        this.m00 = x * m00;
        this.m01 = y * m01;
        this.m10 = x * m10;
        this.m11 = y * m11;
        return this;
    }

    rotateMatrix(matrix, angle) {
        const m00 = matrix.m00, m01 = matrix.m01;
        const m10 = matrix.m10, m11 = matrix.m11;

        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = c * m00 - s * m01;
        this.m01 = c * m01 + s * m00;
        this.m10 = c * m10 - s * m11;
        this.m11 = c * m11 + s * m10;
        return this;
    }

    transpose() {
        return this.transposeMatrix(this);
    }

    inverse() {
        return this.inverseMatrix(this);
    }

    scale(x, y) {
        return this.scaleMatrix(this, x, y);
    }

    rotate(angle) {
        return this.rotateMatrix(this, angle);
    }
}

export class Matrix3 {
    static identity() {
        return new Matrix3(1, 0, 0,
                           0, 1, 0,
                           0, 0, 1);
    }

    static scale(x, y, z) {
        return new Matrix3(x, 0, 0,
                           0, y, 0,
                           0, 0, z);
    }

    static translate(x, y) {
        return new Matrix3(1, 0, x,
                           0, 1, y,
                           0, 0, 1);
    }

    static rotateX(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix3(1, 0, 0,
                           0, c, -s,
                           0, s, c);
    }

    static rotateY(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix3(c, 0, s,
                           0, 1, 0,
                          -s, 0, c);
    }

    static rotateZ(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix3(c, -s, 0,
                           s,  c, 0,
                           0,  0, 1);
    }

    static rotate(axis_x, axis_y, axis_z, angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        const length = Math.sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
        
        if (length == 0) {
            throw new Error("axis vector length must be non zero!");
        }

        const factor = 1 / length;
        const x = factor * axis_x;
        const y = factor * axis_y;
        const z = factor * axis_z;

        const xy = x*y, xs = x*s, xx = x*x;
        const xz = x*z, ys = y*s, yy = y*y;
        const yz = y*z, zs = z*s, zz = z*z;

        const b00 = c + xx*one_minus_c;
        const b10 = xy*one_minus_c - zs;
        const b20 = xz*one_minus_c + ys;

        const b01 = xy*one_minus_c + zs;
        const b11 = yy*one_minus_c + c;
        const b21 = yz*one_minus_c - xs;

        const b02 = xz*one_minus_c - ys;
        const b12 = yz*one_minus_c + xs;
        const b22 = zz*one_minus_c + c;

        return new Matrix3(b00, b10, b20,
                           b01, b11, b21,
                           b02, b12, b22);
    }

    constructor(
        m00=1, m10=0, m20=0,
        m01=0, m11=1, m21=0,
        m02=0, m12=0, m22=1
    ) {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;

        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;

        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
    }

    clone() {
        return new Matrix3(this.m00, this.m10, this.m20,
                           this.m01, this.m11, this.m21,
                           this.m02, this.m12, this.m22);
    }

    determinant() {
        return this.m00 * (this.m11 * this.m22 - this.m12 * this.m21)
             - this.m01 * (this.m10 * this.m22 - this.m12 * this.m20)
             + this.m02 * (this.m10 * this.m21 - this.m11 * this.m20);
    }

    isIdentity() {
        return this.m00 === 1 && this.m01 === 0 && this.m02 === 0
        &&     this.m10 === 0 && this.m11 === 1 && this.m12 === 0
        &&     this.m20 === 0 && this.m21 === 0 && this.m22 === 1;
    }

    equals(other) {
        return Math.abs(this.m00 - other.m00) <= EPSILON && Math.abs(this.m10 - other.m10) <= EPSILON
        &&     Math.abs(this.m20 - other.m20) <= EPSILON && Math.abs(this.m01 - other.m01) <= EPSILON
        &&     Math.abs(this.m11 - other.m11) <= EPSILON && Math.abs(this.m21 - other.m21) <= EPSILON
        &&     Math.abs(this.m02 - other.m02) <= EPSILON && Math.abs(this.m12 - other.m12) <= EPSILON
        &&     Math.abs(this.m22 - other.m22) <= EPSILON;
    }

    identity() {
        this.m00 = 1;
        this.m01 = 0;
        this.m02 = 0;
        this.m10 = 0;
        this.m11 = 1;
        this.m12 = 0;
        this.m20 = 0;
        this.m21 = 0;
        this.m22 = 1;
        return this;
    }

    set(m00, m01, m02, m10, m11, m12, m20, m21, m22) {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
        return this;
    }

    copy(other) {
        this.m00 = other.m00;
        this.m01 = other.m01;
        this.m02 = other.m02;
        this.m10 = other.m10;
        this.m11 = other.m11;
        this.m12 = other.m12;
        this.m20 = other.m20;
        this.m21 = other.m21;
        this.m22 = other.m22;
        return this;
    }

    transposeMatrix(matrix) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22;

        this.m00 = m00;
        this.m01 = m10;
        this.m02 = m20;
        this.m10 = m01;
        this.m11 = m11;
        this.m12 = m21;
        this.m20 = m02;
        this.m21 = m12;
        this.m22 = m22;
        return this;
    }

    mulMatrices(matrix0, matrix1) {
        const a00 = matrix0.m00, a01 = matrix0.m01, a02 = matrix0.m02;
        const a10 = matrix0.m10, a11 = matrix0.m11, a12 = matrix0.m12;
        const a20 = matrix0.m20, a21 = matrix0.m21, a22 = matrix0.m22;

        const b00 = matrix1.m00, b01 = matrix1.m01, b02 = matrix1.m02;
        const b10 = matrix1.m10, b11 = matrix1.m11, b12 = matrix1.m12;
        const b20 = matrix1.m20, b21 = matrix1.m21, b22 = matrix1.m22;

        this.m00 = a00 * b00 + a01 * b10 + a02 * b20;
        this.m01 = a00 * b01 + a01 * b11 + a02 * b21;
        this.m02 = a00 * b02 + a01 * b12 + a02 * b22;
        this.m10 = a10 * b00 + a11 * b10 + a12 * b20;
        this.m11 = a10 * b01 + a11 * b11 + a12 * b21;
        this.m12 = a10 * b02 + a11 * b12 + a12 * b22;
        this.m20 = a20 * b00 + a21 * b10 + a22 * b20;
        this.m21 = a20 * b01 + a21 * b11 + a22 * b21;
        this.m22 = a20 * b02 + a21 * b12 + a22 * b22;
        return this;
    }

    inverseMatrix(matrix) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22;

        const cofactor00 = m11 * m22 - m12 * m21;
        const cofactor10 = m02 * m21 - m01 * m22;
        const cofactor20 = m01 * m12 - m02 * m11;

        const cofactor01 = m12 * m20 - m10 * m22;
        const cofactor11 = m00 * m22 - m02 * m20;
        const cofactor21 = m02 * m10 - m00 * m12;

        const cofactor02 = m10 * m21 - m11 * m20;
        const cofactor12 = m01 * m20 - m00 * m21;
        const cofactor22 = m00 * m11 - m01 * m10;

        const determinant = m00 * cofactor00 + m10 * cofactor10 + m20 * cofactor20;

        if (determinant === 0) {
            throw new Error("Determinant of the matrix cannot be zero!");
        }

        const factor = 1 / determinant;

        this.m00 = cofactor00 * factor;
        this.m01 = cofactor10 * factor;
        this.m02 = cofactor20 * factor;

        this.m10 = cofactor01 * factor;
        this.m11 = cofactor11 * factor;
        this.m12 = cofactor21 * factor;

        this.m20 = cofactor02 * factor;
        this.m21 = cofactor12 * factor;
        this.m22 = cofactor22 * factor;

        return this;
    }

    scaleMatrix(matrix, x, y, z) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22;

        this.m00 = x * m00;
        this.m01 = y * m01;
        this.m02 = z * m02;
        this.m10 = x * m10;
        this.m11 = y * m11;
        this.m12 = z * m12;
        this.m20 = x * m20;
        this.m21 = y * m21;
        this.m22 = z * m22;
        return this;
    }

    translateMatrix(matrix, x, y, z) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22;

        this.m03 = x * m00 + y * m01 + z * m02;
        this.m13 = x * m10 + y * m11 + z * m12;
        this.m23 = x * m20 + y * m21 + z * m22;
        return this;
    }

    rotateMatrixX(matrix, angle) {
        const m01 = matrix.m01, m02 = matrix.m02;
        const m11 = matrix.m11, m12 = matrix.m12;
        const m21 = matrix.m21, m22 = matrix.m22;

        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m01 = c * m01 - s * m02;
        this.m02 = c * m02 + s * m01;
        this.m11 = c * m11 - s * m12;
        this.m12 = c * m12 + s * m11;
        this.m21 = c * m21 - s * m22;
        this.m22 = c * m22 + s * m21;
        return this;
    }

    rotateMatrixY(matrix, angle) {
        const m00 = matrix.m00, m02 = matrix.m02;
        const m10 = matrix.m10, m12 = matrix.m12;
        const m20 = matrix.m20, m22 = matrix.m22;

        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = c * m00 + s * m02;
        this.m02 = c * m02 - s * m00;
        this.m10 = c * m10 + s * m12;
        this.m12 = c * m12 - s * m10;
        this.m20 = c * m20 + s * m22;
        this.m22 = c * m22 - s * m20;
        return this;
    }

    rotateMatrixZ(matrix, angle) {
        const m00 = matrix.m00, m01 = matrix.m01;
        const m10 = matrix.m10, m11 = matrix.m11;
        const m20 = matrix.m20, m21 = matrix.m21;

        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = c * m00 - s * m01;
        this.m01 = c * m01 + s * m00;
        this.m10 = c * m10 - s * m11;
        this.m11 = c * m11 + s * m10;
        this.m20 = c * m20 - s * m21;
        this.m21 = c * m21 + s * m20;
        return this;
    }

    rotateMatrix(matrix, axis_x, axis_y, axis_z, angle) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22;

        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        const length = Math.sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);

        if (length == 0) {
            throw new Error("axis vector length must be non zero!");
        }

        const factor = 1 / length;
        const x = factor * axis_x;
        const y = factor * axis_y;
        const z = factor * axis_z;

        const xy = x * y, xs = x * s, xx = x * x;
        const xz = x * z, ys = y * s, yy = y * y;
        const yz = y * z, zs = z * s, zz = z * z;

        const b00 = c + xx * one_minus_c;
        const b10 = xy * one_minus_c - zs;
        const b20 = xz * one_minus_c + ys;

        const b01 = xy * one_minus_c + zs;
        const b11 = yy * one_minus_c + c;
        const b21 = yz * one_minus_c - xs;

        const b02 = xz * one_minus_c - ys;
        const b12 = yz * one_minus_c + xs;
        const b22 = zz * one_minus_c + c;

        this.m00 = m10 * b10 + m20 * b20 + m00 * b00;
        this.m01 = m00 * b01 + m20 * b21 + m10 * b11;
        this.m02 = m00 * b02 + m10 * b12 + m20 * b22;
        this.m10 = m11 * b10 + m21 * b20 + m01 * b00;
        this.m11 = m01 * b01 + m21 * b21 + m11 * b11;
        this.m12 = m01 * b02 + m11 * b12 + m21 * b22;
        this.m20 = m12 * b10 + m22 * b20 + m02 * b00;
        this.m21 = m02 * b01 + m22 * b21 + m12 * b11;
        this.m22 = m02 * b02 + m12 * b12 + m22 * b22;
        return this;
    }
    
    transpose() {
        return this.transposeMatrix(this);
    }

    inverse() {
        return this.inverseMatrix(this);
    }

    scale(x, y, z) {
        return this.scaleMatrix(this, x, y, z);
    }

    translate(x, y, z) {
        return this.translateMatrix(this, x, y, z);
    }

    rotateX(angle) {
        return this.rotateMatrixX(this, angle);
    }

    rotateY(angle) {
        return this.rotateMatrixY(this, angle);
    }

    rotateZ(angle) {
        return this.rotateMatrixZ(this, angle);
    }

    rotate(axis_x, axis_y, axis_z, angle) {
        return this.rotateMatrix(this, axis_x, axis_y, axis_z, angle);
    }
}

export class Matrix4 {
    static identity() {
        return new Matrix4(1, 0, 0, 0,
                           0, 1, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 1);
    }

    static scale(x, y, z) {
        return new Matrix4(x, 0, 0, 0,
                           0, y, 0, 0,
                           0, 0, z, 0,
                           0, 0, 0, 1);
    }

    static translate(x, y, z) {
        return new Matrix4(0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           x, y, z, 1);
    }

    static rotateX(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(1, 0, 0, 0,
                           0, c,-s, 0,
                           0, s, c, 0,
                           0, 0, 0, 1);
    }

    static rotateY(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(c, 0, s, 0,
                           0, 1, 0, 0,
                          -s, 0, c, 0,
                           0, 0, 0, 1);
    }

    static rotateZ(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(c,-s, 0, 0,
                           s, c, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 1);
    }

    static rotate(axis_x, axis_y, axis_z, angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        const length = Math.sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
        
        if (length == 0) {
            throw new Error("axis vector length must be non zero!");
        }

        const factor = 1 / length;
        const x = factor * axis_x;
        const y = factor * axis_y;
        const z = factor * axis_z;

        const xy = x*y, xs = x*s, xx = x*x;
        const xz = x*z, ys = y*s, yy = y*y;
        const yz = y*z, zs = z*s, zz = z*z;

        const b00 = c + xx*one_minus_c;
        const b10 = xy*one_minus_c - zs;
        const b20 = xz*one_minus_c + ys;

        const b01 = xy*one_minus_c + zs;
        const b11 = yy*one_minus_c + c;
        const b21 = yz*one_minus_c - xs;

        const b02 = xz*one_minus_c - ys;
        const b12 = yz*one_minus_c + xs;
        const b22 = zz*one_minus_c + c;

        return new Matrix4(b00, b10, b20, 0,
                           b01, b11, b21, 0,
                           b02, b12, b22, 0,
                           0,   0,   0,   1);
    }

    constructor(
        m00=1, m10=0, m20=0, m30=0,
        m01=0, m11=1, m21=0, m31=0,
        m02=0, m12=0, m22=1, m32=0,
        m03=0, m13=0, m23=0, m33=1,
    ) {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m03 = m03;

        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m13 = m13;

        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
        this.m23 = m23;

        this.m30 = m30;
        this.m31 = m31;
        this.m32 = m32;
        this.m33 = m33;
    }

    clone() {
        return new Matrix4(this.m00, this.m10, this.m20, this.m30,
                           this.m01, this.m11, this.m21, this.m31,
                           this.m02, this.m12, this.m22, this.m32,
                           this.m03, this.m13, this.m23, this.m33);
    }

    determinant() {
        const subfactor0 = this.m22 * this.m33 - this.m23 * this.m32;
        const subfactor1 = this.m21 * this.m33 - this.m23 * this.m31;
        const subfactor2 = this.m21 * this.m32 - this.m22 * this.m31;
        const subfactor3 = this.m20 * this.m33 - this.m23 * this.m30;
        const subfactor4 = this.m20 * this.m32 - this.m22 * this.m30;
        const subfactor5 = this.m20 * this.m31 - this.m21 * this.m30;

        const cofactor00 = this.m11 * subfactor0 - this.m12 * subfactor1 + this.m13 * subfactor2;
        const cofactor01 = this.m10 * subfactor0 - this.m12 * subfactor3 + this.m13 * subfactor4;
        const cofactor02 = this.m10 * subfactor1 - this.m11 * subfactor3 + this.m13 * subfactor5;
        const cofactor03 = this.m10 * subfactor2 - this.m11 * subfactor4 + this.m12 * subfactor5;

        return this.m00 * cofactor00 - this.m01 * cofactor01
             + this.m02 * cofactor02 - this.m03 * cofactor03;
    }

    isIdentity() {
        return this.m00 === 1 && this.m01 === 0 && this.m02 === 0 && this.m03 === 0
        &&     this.m10 === 0 && this.m11 === 1 && this.m12 === 0 && this.m13 === 0
        &&     this.m20 === 0 && this.m21 === 0 && this.m22 === 1 && this.m23 === 0
        &&     this.m30 === 0 && this.m31 === 0 && this.m32 === 0 && this.m33 === 1;
    }

    equals(other) {
        return Math.abs(this.m00 - other.m00) <= EPSILON && Math.abs(this.m10 - other.m10) <= EPSILON
        &&     Math.abs(this.m20 - other.m20) <= EPSILON && Math.abs(this.m30 - other.m30) <= EPSILON
        &&     Math.abs(this.m01 - other.m01) <= EPSILON && Math.abs(this.m11 - other.m11) <= EPSILON
        &&     Math.abs(this.m21 - other.m21) <= EPSILON && Math.abs(this.m31 - other.m31) <= EPSILON
        &&     Math.abs(this.m02 - other.m02) <= EPSILON && Math.abs(this.m12 - other.m12) <= EPSILON
        &&     Math.abs(this.m22 - other.m22) <= EPSILON && Math.abs(this.m32 - other.m32) <= EPSILON
        &&     Math.abs(this.m03 - other.m03) <= EPSILON && Math.abs(this.m13 - other.m13) <= EPSILON
        &&     Math.abs(this.m23 - other.m23) <= EPSILON && Math.abs(this.m33 - other.m33) <= EPSILON;
    }
    
    identity() {
        this.m00 = 1;
        this.m01 = 0;
        this.m02 = 0;
        this.m03 = 0;
        this.m10 = 0;
        this.m11 = 1;
        this.m12 = 0;
        this.m13 = 0;
        this.m20 = 0;
        this.m21 = 0;
        this.m22 = 1;
        this.m23 = 0;
        this.m30 = 0;
        this.m31 = 0;
        this.m32 = 0;
        this.m33 = 1;
        return this;
    }

    set(m00, m01, m02, m03, m10, m11, m12, m13,
        m20, m21, m22, m23, m30, m31, m32, m33)
    {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m03 = m03;
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m13 = m13;
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
        this.m23 = m23;
        this.m30 = m30;
        this.m31 = m31;
        this.m32 = m32;
        this.m33 = m33;
        return this;
    }

    
    setCols(m00, m01, m02, m03, m10, m11, m12, m13,
            m20, m21, m22, m23, m30, m31, m32, m33)
    {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m03 = m03;
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m13 = m13;
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
        this.m23 = m23;
        this.m30 = m30;
        this.m31 = m31;
        this.m32 = m32;
        this.m33 = m33;
        return this;
    }

    
    setRows(m00, m10, m20, m30, m01, m11, m21, m31,
            m02, m12, m22, m32, m03, m13, m23, m33)
    {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m03 = m03;
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m13 = m13;
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
        this.m23 = m23;
        this.m30 = m30;
        this.m31 = m31;
        this.m32 = m32;
        this.m33 = m33;
        return this;
    }


    copy(other) {
        this.m00 = other.m00;
        this.m01 = other.m01;
        this.m02 = other.m02;
        this.m03 = other.m03;
        this.m10 = other.m10;
        this.m11 = other.m11;
        this.m12 = other.m12;
        this.m13 = other.m13;
        this.m20 = other.m20;
        this.m21 = other.m21;
        this.m22 = other.m22;
        this.m23 = other.m23;
        this.m30 = other.m30;
        this.m31 = other.m31;
        this.m32 = other.m32;
        this.m33 = other.m33;
        return this;
    }

    transposeMatrix(matrix) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02, m03 = matrix.m03;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12, m13 = matrix.m13;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22, m23 = matrix.m23;
        const m30 = matrix.m30, m31 = matrix.m31, m32 = matrix.m32, m33 = matrix.m33;

        this.m00 = m00;
        this.m01 = m10;
        this.m02 = m20;
        this.m03 = m30;
        this.m10 = m01;
        this.m11 = m11;
        this.m12 = m21;
        this.m13 = m31;
        this.m20 = m02;
        this.m21 = m12;
        this.m22 = m22;
        this.m23 = m32;
        this.m30 = m03;
        this.m31 = m13;
        this.m32 = m23;
        this.m33 = m33;
        return this;
    }

    mulMatrices(matrix0, matrix1) {
        const a00 = matrix0.m00, a01 = matrix0.m01, a02 = matrix0.m02, a03 = matrix0.m03;
        const a10 = matrix0.m10, a11 = matrix0.m11, a12 = matrix0.m12, a13 = matrix0.m13;
        const a20 = matrix0.m20, a21 = matrix0.m21, a22 = matrix0.m22, a23 = matrix0.m23;
        const a30 = matrix0.m30, a31 = matrix0.m31, a32 = matrix0.m32, a33 = matrix0.m33;

        const b00 = matrix1.m00, b01 = matrix1.m01, b02 = matrix1.m02, b03 = matrix1.m03;
        const b10 = matrix1.m10, b11 = matrix1.m11, b12 = matrix1.m12, b13 = matrix1.m13;
        const b20 = matrix1.m20, b21 = matrix1.m21, b22 = matrix1.m22, b23 = matrix1.m23;
        const b30 = matrix1.m30, b31 = matrix1.m31, b32 = matrix1.m32, b33 = matrix1.m33;

        this.m00 = a00 * b00 + a10 * b01 + a20 * b02 + a30 * b03;
        this.m01 = a00 * b10 + a10 * b11 + a20 * b12 + a30 * b13;
        this.m02 = a00 * b20 + a10 * b21 + a20 * b22 + a30 * b23;
        this.m03 = a00 * b30 + a10 * b31 + a20 * b32 + a30 * b33;
        this.m10 = a01 * b00 + a11 * b01 + a21 * b02 + a31 * b03;
        this.m11 = a01 * b10 + a11 * b11 + a21 * b12 + a31 * b13;
        this.m12 = a01 * b20 + a11 * b21 + a21 * b22 + a31 * b23;
        this.m13 = a01 * b30 + a11 * b31 + a21 * b32 + a31 * b33;
        this.m20 = a02 * b00 + a12 * b01 + a22 * b02 + a32 * b03;
        this.m21 = a02 * b10 + a12 * b11 + a22 * b12 + a32 * b13;
        this.m22 = a02 * b20 + a12 * b21 + a22 * b22 + a32 * b23;
        this.m23 = a02 * b30 + a12 * b31 + a22 * b32 + a32 * b33;
        this.m30 = a03 * b00 + a13 * b01 + a23 * b02 + a33 * b03;
        this.m31 = a03 * b10 + a13 * b11 + a23 * b12 + a33 * b13;
        this.m32 = a03 * b20 + a13 * b21 + a23 * b22 + a33 * b23;
        this.m33 = a03 * b30 + a13 * b31 + a23 * b32 + a33 * b33;
        return this;
    }

    inverseMatrix(matrix) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02, m03 = matrix.m03;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12, m13 = matrix.m13;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22, m23 = matrix.m23;
        const m30 = matrix.m30, m31 = matrix.m31, m32 = matrix.m32, m33 = matrix.m33;

        const subfactor0  = m22 * m33 - m23 * m32;
        const subfactor1  = m12 * m33 - m23 * m31;
        const subfactor2  = m12 * m32 - m22 * m31;
        const subfactor3  = m02 * m33 - m23 * m30;
        const subfactor4  = m02 * m32 - m22 * m30;
        const subfactor5  = m02 * m31 - m21 * m30;

        const subfactor6  = m20 * m13 - m03 * m12;
        const subfactor7  = m10 * m13 - m03 * m11;
        const subfactor8  = m10 * m12 - m02 * m11;
        const subfactor9  = m00 * m13 - m03 * m10;
        const subfactor10 = m00 * m12 - m02 * m10;
        const subfactor11 = m00 * m11 - m01 * m10;

        const cofactor00 = m11 * subfactor0  - m12 * subfactor1  + m13 * subfactor2;
        const cofactor01 = m21 * subfactor3  - m13 * subfactor4  - m10 * subfactor0;
        const cofactor02 = m01 * subfactor1  - m11 * subfactor3  + m13 * subfactor5;
        const cofactor03 = m11 * subfactor4  - m12 * subfactor5  - m10 * subfactor2;

        const cofactor10 = m20 * subfactor1  - m03 * subfactor2  - m01 * subfactor0;
        const cofactor11 = m00 * subfactor0  - m02 * subfactor3  + m03 * subfactor4;
        const cofactor12 = m10 * subfactor3  - m03 * subfactor5  - m00 * subfactor1;
        const cofactor13 = m00 * subfactor2  - m01 * subfactor4  + m02 * subfactor5;

        const cofactor20 = m13 * subfactor6  - m32 * subfactor7  + m33 * subfactor8;
        const cofactor21 = m23 * subfactor9  - m33 * subfactor10 - m30 * subfactor6;
        const cofactor22 = m03 * subfactor7  - m31 * subfactor9  + m33 * subfactor11;
        const cofactor23 = m13 * subfactor10 - m32 * subfactor11 - m30 * subfactor8;

        const cofactor30 = m22 * subfactor7  - m23 * subfactor8  - m21 * subfactor6;
        const cofactor31 = m02 * subfactor6  - m22 * subfactor9  + m23 * subfactor10;
        const cofactor32 = m12 * subfactor9  - m23 * subfactor11 - m20 * subfactor7;
        const cofactor33 = m02 * subfactor8  - m21 * subfactor10 + m22 * subfactor11;

        const determinant = m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02 + m03 * cofactor03;

        if (determinant === 0) {
            throw new Error("Determinant of the cofactor matrix cannot be zero!");
        }

        const factor = 1 / determinant;

        this.m00 = cofactor00 * factor;
        this.m01 = cofactor10 * factor;
        this.m02 = cofactor20 * factor;
        this.m03 = cofactor30 * factor;

        this.m10 = cofactor01 * factor;
        this.m11 = cofactor11 * factor;
        this.m12 = cofactor21 * factor;
        this.m13 = cofactor31 * factor;

        this.m20 = cofactor02 * factor;
        this.m21 = cofactor12 * factor;
        this.m22 = cofactor22 * factor;
        this.m23 = cofactor32 * factor;

        this.m30 = cofactor03 * factor;
        this.m31 = cofactor13 * factor;
        this.m32 = cofactor23 * factor;
        this.m33 = cofactor33 * factor;

        return this;
    }

    scaleMatrix(matrix, x, y, z) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22;
        const m30 = matrix.m30, m31 = matrix.m31, m32 = matrix.m32;

        this.m00 = x * m00;
        this.m01 = y * m01;
        this.m02 = z * m02;
        this.m10 = x * m10;
        this.m11 = y * m11;
        this.m12 = z * m12;
        this.m20 = x * m20;
        this.m21 = y * m21;
        this.m22 = z * m22;
        this.m30 = x * m30;
        this.m31 = y * m31;
        this.m32 = z * m32;
        return this;
    }

    translateMatrix(matrix, x, y, z) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02, m03 = matrix.m03;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12, m13 = matrix.m13;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22, m23 = matrix.m23;
        const m30 = matrix.m30, m31 = matrix.m31, m32 = matrix.m32, m33 = matrix.m33;

        this.m00 = m00 + x * m03;
        this.m10 = m10 + x * m13;
        this.m20 = m20 + x * m23;
        this.m30 = m30 + x * m33;

        this.m01 = m01 + y * m03;
        this.m11 = m11 + y * m13;
        this.m21 = m21 + y * m23;
        this.m31 = m31 + y * m33;
        
        this.m02 = m02 + z * m03;
        this.m12 = m12 + z * m13;
        this.m22 = m22 + z * m23;
        this.m32 = m32 + z * m33;

        this.m03 = m03;
        this.m13 = m13;
        this.m23 = m23;
        this.m33 = m33;

        // this.m00 = x * m03 + m00;
        // this.m10 = x * m13 + m10;
        // this.m20 = x * m23 + m20;
        // this.m30 = x * m33 + m30;
        // this.m01 = y * m03 + m01;
        // this.m11 = y * m13 + m11;
        // this.m21 = y * m23 + m21;
        // this.m31 = y * m33 + m31;
        // this.m02 = z * m03 + m02;
        // this.m12 = z * m13 + m12;
        // this.m22 = z * m23 + m22;
        // this.m32 = z * m03 + m32;
        // this.m03 = m03;
        // this.m13 = m13;
        // this.m23 = m23;
        // this.m33 = m33;

        // this.m03 = x * m00 + y * m01 + z * m02 + m03;
        // this.m13 = x * m10 + y * m11 + z * m12 + m13;
        // this.m23 = x * m20 + y * m21 + z * m22 + m23;
        // this.m33 = x * m30 + y * m31 + z * m32 + m33;
        
        // this.m30 = x * m00 + y * m01 + z * m02 + m03;
        // this.m31 = x * m10 + y * m11 + z * m12 + m13;
        // this.m32 = x * m20 + y * m21 + z * m22 + m23;
        // this.m33 = x * m30 + y * m31 + z * m32 + m33;
        return this;
    }

    rotateMatrixX(matrix, angle) {
        const m01 = matrix.m01, m02 = matrix.m02;
        const m11 = matrix.m11, m12 = matrix.m12;
        const m21 = matrix.m21, m22 = matrix.m22;
        const m31 = matrix.m31, m32 = matrix.m32;

        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m01 = c * m01 - s * m02;
        this.m02 = c * m02 + s * m01;
        this.m11 = c * m11 - s * m12;
        this.m12 = c * m12 + s * m11;
        this.m21 = c * m21 - s * m22;
        this.m22 = c * m22 + s * m21;
        this.m31 = c * m31 - s * m32;
        this.m32 = c * m32 + s * m31;
        return this;
    }

    rotateMatrixY(matrix, angle) {
        const m00 = matrix.m00, m02 = matrix.m02;
        const m10 = matrix.m10, m12 = matrix.m12;
        const m20 = matrix.m20, m22 = matrix.m22;
        const m30 = matrix.m30, m32 = matrix.m32;

        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = c * m00 + s * m02;
        this.m02 = c * m02 - s * m00;
        this.m10 = c * m10 + s * m12;
        this.m12 = c * m12 - s * m10;
        this.m20 = c * m20 + s * m22;
        this.m22 = c * m22 - s * m20;
        this.m30 = c * m30 + s * m32;
        this.m32 = c * m32 - s * m30;
        return this;
    }

    rotateMatrixZ(matrix, angle) {
        const m00 = matrix.m00, m01 = matrix.m01;
        const m10 = matrix.m10, m11 = matrix.m11;
        const m20 = matrix.m20, m21 = matrix.m21;
        const m30 = matrix.m30, m31 = matrix.m31;

        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = c * m00 - s * m01;
        this.m01 = c * m01 + s * m00;
        this.m10 = c * m10 - s * m11;
        this.m11 = c * m11 + s * m10;
        this.m20 = c * m20 - s * m21;
        this.m21 = c * m21 + s * m20;
        this.m30 = c * m30 - s * m31;
        this.m31 = c * m31 + s * m30;
        return this;
    }

    rotateMatrix(matrix, axis_x, axis_y, axis_z, angle) {
        const m00 = matrix.m00, m01 = matrix.m01, m02 = matrix.m02, m03 = matrix.m03;
        const m10 = matrix.m10, m11 = matrix.m11, m12 = matrix.m12, m13 = matrix.m13;
        const m20 = matrix.m20, m21 = matrix.m21, m22 = matrix.m22, m23 = matrix.m23;

        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        const length = Math.sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);

        if (length == 0) {
            throw new Error("axis vector length must be non zero!");
        }

        const factor = 1 / length;
        const x = factor * axis_x;
        const y = factor * axis_y;
        const z = factor * axis_z;

        const xy = x * y, xs = x * s, xx = x * x;
        const xz = x * z, ys = y * s, yy = y * y;
        const yz = y * z, zs = z * s, zz = z * z;

        const b00 = c + xx * one_minus_c;
        const b10 = xy * one_minus_c - zs;
        const b20 = xz * one_minus_c + ys;

        const b01 = xy * one_minus_c + zs;
        const b11 = yy * one_minus_c + c;
        const b21 = yz * one_minus_c - xs;

        const b02 = xz * one_minus_c - ys;
        const b12 = yz * one_minus_c + xs;
        const b22 = zz * one_minus_c + c;

        this.m00 = m10 * b10 + m20 * b20 + m00 * b00;
        this.m01 = m00 * b01 + m20 * b21 + m10 * b11;
        this.m02 = m00 * b02 + m10 * b12 + m20 * b22;
        this.m10 = m11 * b10 + m21 * b20 + m01 * b00;
        this.m11 = m01 * b01 + m21 * b21 + m11 * b11;
        this.m12 = m01 * b02 + m11 * b12 + m21 * b22;
        this.m20 = m12 * b10 + m22 * b20 + m02 * b00;
        this.m21 = m02 * b01 + m22 * b21 + m12 * b11;
        this.m22 = m02 * b02 + m12 * b12 + m22 * b22;
        this.m30 = m13 * b10 + m23 * b20 + m03 * b00;
        this.m31 = m03 * b01 + m23 * b21 + m13 * b11;
        this.m32 = m03 * b02 + m13 * b12 + m23 * b22;
        return this;
    }

    transpose() {
        return this.transposeMatrix(this);
    }

    inverse() {
        return this.inverseMatrix(this);
    }

    scale(x, y, z) {
        return this.scaleMatrix(this, x, y, z);
    }

    translate(x, y, z) {
        return this.translateMatrix(this, x, y, z);
    }

    rotateX(angle) {
        return this.rotateMatrixX(this, angle);
    }

    rotateY(angle) {
        return this.rotateMatrixY(this, angle);
    }

    rotateZ(angle) {
        return this.rotateMatrixZ(this, angle);
    }

    rotate(axis_x, axis_y, axis_z, angle) {
        return this.rotateMatrix(this, axis_x, axis_y, axis_z, angle);
    }
}

export function clamp(value, min, max) {
    if (value > max) return max;
    if (value < min) return min;
    return value;
}

export function gcf(a, b) {
    if (a < b) {
        return gcf(b, a);
    }

    if (Math.abs(b) < 0.001) {
        return a;
    }

    return gcf(b, a - Math.floor(a / b) * b);
}

export function lcm(a, b) {
    return a * b / gcf(a, b);
}

export function lerp(a, b, t) {
    return (b - a) * t + a;
}

export function smallestGridSize(size_x, size_y, size_z) {
    return gcf(size_x, gcf(size_y, size_z));
}

export function closestEvenlyDivisibleGridSize(size_x, size_y, size_z, grid_size) {
    if (size_x == size_y && size_x == size_z) {
        const outside_volume = size_x / grid_size;
        const inside_volume = Math.floor(outside_volume);
        return (inside_volume / outside_volume) * grid_size;
    }

    const smallest_grid_size = gcf(size_x, gcf(size_y, size_z));

    if (grid_size < smallest_grid_size) {
        return Math.min(Math.floor(smallest_grid_size / grid_size), 1) * grid_size;
    }

    return Math.max(Math.floor(grid_size / smallest_grid_size), 1) * smallest_grid_size;
}

export function bestDivisibleGridSize(size_x, size_y, size_z, division_count) {
    const grid_size = Math.max(size_x, size_y, size_z) / division_count;
    return closestEvenlyDivisibleGridSize(size_x, size_y, size_z, grid_size);
}

export function distanceFromOrigin(x, y, z) {
    return Math.sqrt(x*x + y*y + z*z);
}

export function distanceSquared(x0, y0, z0, x1, y1, z1) {
    const dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
    return dx * dx + dy * dy + dz * dz;
}

export class Vertex {
    constructor(
        position = new Vector3(0, 0, 0),
        normal = new Vector3(0, 0, 0),
        colour = new Vector4(0, 0, 0, 1),
        texture_coords = new Vector2(0, 0),
        texture_id = 0
    ) {
        this.position = position;
        this.normal = normal;
        this.colour = colour;
        this.texture_coords = texture_coords;
        this.texture_id = texture_id;
    }
}

export class TriangleMesh {
    constructor(vertices = [], indices = []) {
        this.vertices = vertices;
        this.indices = indices;
    }
}

export class PolygonIndices {
    constructor(indices = [], plane = new Plane(), colour = null) {
        this.indices = indices;
        this.plane = plane;
        this.colour = colour;
    }
}

export class PolygonMesh {
    constructor(vertices = [], faces = []) {
        this.vertices = vertices;
        this.faces = faces;
    }
}

function getLowestSubdivision(node0, node1, node2, node3) {
    let lowest = node0.subdivision_count < node1.subdivision_count ? node0 : node1;
    lowest = lowest.subdivision_count < node2.subdivision_count ? lowest : node2;
    lowest = lowest.subdivision_count < node3.subdivision_count ? lowest : node3;
    return lowest;
}

function interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, x, y, z) {
    const dx = (x - min_x) / (max_x - min_x), one_minus_dx = 1 - dx;
    const dy = (y - min_y) / (max_y - min_y), one_minus_dy = 1 - dy;
    const dz = (z - min_z) / (max_z - min_z), one_minus_dz = 1 - dz;

    const c000 = field.calculateSignedDistance(min_x, min_y, min_z);
    const c100 = field.calculateSignedDistance(max_x, min_y, min_z);
    const c110 = field.calculateSignedDistance(max_x, max_y, min_z);
    const c010 = field.calculateSignedDistance(min_x, max_y, min_z);
    const c001 = field.calculateSignedDistance(min_x, min_y, max_z);
    const c101 = field.calculateSignedDistance(max_x, min_y, max_z);
    const c111 = field.calculateSignedDistance(max_x, max_y, max_z);
    const c011 = field.calculateSignedDistance(min_x, max_y, max_z);

    const x00 = c000 * one_minus_dx + c100 * dx;
    const x01 = c001 * one_minus_dx + c101 * dx;
    const x10 = c010 * one_minus_dx + c110 * dx;
    const x11 = c011 * one_minus_dx + c111 * dx;

    const y0 = x00 * one_minus_dy + x10 * dy;
    const y1 = x01 * one_minus_dy + x11 * dy;

    return y0 * one_minus_dz + y1 * dz;
}

class SampleEdge {
    constructor(x0, y0, z0, sample0, x1, y1, z1, sample1) {
        this.x0 = x0;
        this.y0 = y0;
        this.z0 = z0;
        this.x1 = x1;
        this.y1 = y1;
        this.z1 = z1;
        this.sample0 = sample0;
        this.sample1 = sample1;
    }

    set(x0, y0, z0, sample0, x1, y1, z1, sample1) {
        this.x0 = x0;
        this.y0 = y0;
        this.z0 = z0;
        this.x1 = x1;
        this.y1 = y1;
        this.z1 = z1;
        this.sample0 = sample0;
        this.sample1 = sample1;
    }
}

export class SampleOctreeRoot {
    constructor(min_x, min_y, min_z, max_x, max_y, max_z, subdivision_count = 1) {
        this.kind  = SAMPLE_OCTREE_KIND_ROOT;
        this.min_x = min_x;
        this.min_y = min_y;
        this.min_z = min_z;
        this.max_x = max_x;
        this.max_y = max_y;
        this.max_z = max_z;
        this.subdivision_count = subdivision_count;
        this.node000 = null;
        this.node100 = null;
        this.node110 = null;
        this.node010 = null;
        this.node001 = null;
        this.node101 = null;
        this.node111 = null;
        this.node011 = null;
        this.sample000 = null;
        this.sample100 = null;
        this.sample110 = null;
        this.sample010 = null;
        this.sample001 = null;
        this.sample101 = null;
        this.sample111 = null;
        this.sample011 = null;
    }

    sampleSignedDistanceField(field, merge_threshold = null) {
        const min_x = this.min_x, max_x = this.max_x, center_x = (max_x + min_x) / 2;
        const min_y = this.min_y, max_y = this.max_y, center_y = (max_y + min_y) / 2;
        const min_z = this.min_z, max_z = this.max_z, center_z = (max_z + min_z) / 2;

        const center_sample = field.calculateSignedDistance(center_x, center_y, center_z);
        const center_right_sample = field.calculateSignedDistance(max_x, center_y, center_z);
        const center_left_sample = field.calculateSignedDistance(min_x, center_y, center_z);
        const center_near_sample = field.calculateSignedDistance(center_x, min_y, center_z);
        const center_far_sample = field.calculateSignedDistance(center_x, max_y, center_z);
        const center_top_sample = field.calculateSignedDistance(center_x, center_y, max_z);
        const center_bottom_sample = field.calculateSignedDistance(center_x, center_y, min_z);

        if (merge_threshold != null) {
            const center_interp = interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, center_x, center_y, center_z);
            const center_right_interp = interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, max_x, center_y, center_z);
            const center_left_interp = interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, min_x, center_y, center_z);
            const center_near_interp = interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, center_x, min_y, center_z);
            const center_far_interp = interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, center_x, max_y, center_z);
            const center_top_interp = interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, center_x, center_y, max_z);
            const center_bottom_interp = interpolatePoint(field, min_x, min_y, min_z, max_x, max_y, max_z, center_x, center_y, min_z);
            
            const accurate_bilinear_interpolation = Math.abs(center_sample - center_interp) <= merge_threshold
            &&                                    Math.abs(center_right_sample - center_right_interp) <= merge_threshold
            &&                                    Math.abs(center_left_sample - center_left_interp) <= merge_threshold
            &&                                    Math.abs(center_near_sample - center_near_interp) <= merge_threshold
            &&                                    Math.abs(center_far_sample - center_far_interp) <= merge_threshold
            &&                                    Math.abs(center_top_sample - center_top_interp) <= merge_threshold
            &&                                    Math.abs(center_bottom_sample - center_bottom_interp) <= merge_threshold;

            if (accurate_bilinear_interpolation) {
                const result = new SampleOctreeLeaf(min_x, min_y, min_z, max_x, max_y, max_z);
                result.sample000 = field.calculateSignedDistance(min_x, min_y, min_z);
                result.sample001 = field.calculateSignedDistance(min_x, min_y, max_z);
                result.sample010 = field.calculateSignedDistance(min_x, max_y, min_z);
                result.sample011 = field.calculateSignedDistance(min_x, max_y, max_z);
                result.sample100 = field.calculateSignedDistance(max_x, min_y, min_z);
                result.sample101 = field.calculateSignedDistance(max_x, min_y, max_z);
                result.sample110 = field.calculateSignedDistance(max_x, max_y, min_z);
                result.sample111 = field.calculateSignedDistance(max_x, max_y, max_z);
                return result.sampleSignedDistanceField(field);
            }
        }

        if (this.subdivision_count == 1) {
            this.node000 = new SampleOctreeLeaf(min_x, min_y, min_z, center_x, center_y, center_z, 0);
            this.node001 = new SampleOctreeLeaf(min_x, min_y, center_z, center_x, center_y, max_z, 0);
            this.node010 = new SampleOctreeLeaf(min_x, center_y, min_z, center_x, max_y, center_z, 0);
            this.node011 = new SampleOctreeLeaf(min_x, center_y, center_z, center_x, max_y, max_z, 0);
            this.node100 = new SampleOctreeLeaf(center_x, min_y, min_z, max_x, center_y, center_z, 0);
            this.node101 = new SampleOctreeLeaf(center_x, min_y, center_z, max_x, center_y, max_z, 0);
            this.node110 = new SampleOctreeLeaf(center_x, center_y, min_z, max_x, max_y, center_z, 0);
            this.node111 = new SampleOctreeLeaf(center_x, center_y, center_z, max_x, max_y, max_z, 0);
        }
        else {
            const subdivision_count = this.subdivision_count - 1;
            this.node000 = new SampleOctreeRoot(min_x, min_y, min_z, center_x, center_y, center_z, subdivision_count);
            this.node001 = new SampleOctreeRoot(min_x, min_y, center_z, center_x, center_y, max_z, subdivision_count);
            this.node010 = new SampleOctreeRoot(min_x, center_y, min_z, center_x, max_y, center_z, subdivision_count);
            this.node011 = new SampleOctreeRoot(min_x, center_y, center_z, center_x, max_y, max_z, subdivision_count);
            this.node100 = new SampleOctreeRoot(center_x, min_y, min_z, max_x, center_y, center_z, subdivision_count);
            this.node101 = new SampleOctreeRoot(center_x, min_y, center_z, max_x, center_y, max_z, subdivision_count);
            this.node110 = new SampleOctreeRoot(center_x, center_y, min_z, max_x, max_y, center_z, subdivision_count);
            this.node111 = new SampleOctreeRoot(center_x, center_y, center_z, max_x, max_y, max_z, subdivision_count);
        }

        const center_min_min_sample = field.calculateSignedDistance(center_x, min_y, min_z);
        const center_max_min_sample = field.calculateSignedDistance(center_x, max_y, min_z);
        const center_max_max_sample = field.calculateSignedDistance(center_x, max_y, max_z);
        const center_min_max_sample = field.calculateSignedDistance(center_x, min_y, max_z);
        const min_center_min_sample = field.calculateSignedDistance(min_x, center_y, min_z);
        const max_center_min_sample = field.calculateSignedDistance(max_x, center_y, min_z);
        const max_center_max_sample = field.calculateSignedDistance(max_x, center_y, max_z);
        const min_center_max_sample = field.calculateSignedDistance(min_x, center_y, max_z);
        const min_min_center_sample = field.calculateSignedDistance(min_x, min_y, center_z);
        const max_min_center_sample = field.calculateSignedDistance(max_x, min_y, center_z);
        const max_max_center_sample = field.calculateSignedDistance(max_x, max_y, center_z);
        const min_max_center_sample = field.calculateSignedDistance(min_x, max_y, center_z);

        // calculate corners
        this.sample000 ??= field.calculateSignedDistance(min_x, min_y, min_z);
        this.sample001 ??= field.calculateSignedDistance(min_x, min_y, max_z);
        this.sample010 ??= field.calculateSignedDistance(min_x, max_y, min_z);
        this.sample011 ??= field.calculateSignedDistance(min_x, max_y, max_z);
        this.sample100 ??= field.calculateSignedDistance(max_x, min_y, min_z);
        this.sample101 ??= field.calculateSignedDistance(max_x, min_y, max_z);
        this.sample110 ??= field.calculateSignedDistance(max_x, max_y, min_z);
        this.sample111 ??= field.calculateSignedDistance(max_x, max_y, max_z);

        // distibute corners
        this.node000.sample000 = this.sample000;
        this.node001.sample001 = this.sample001;
        this.node010.sample010 = this.sample010;
        this.node011.sample011 = this.sample011;
        this.node100.sample100 = this.sample100;
        this.node101.sample101 = this.sample101;
        this.node110.sample110 = this.sample110;
        this.node111.sample111 = this.sample111;

        // distibute center
        this.node000.sample111 = center_sample;
        this.node001.sample110 = center_sample;
        this.node010.sample101 = center_sample;
        this.node011.sample100 = center_sample;
        this.node100.sample011 = center_sample;
        this.node101.sample010 = center_sample;
        this.node110.sample001 = center_sample;
        this.node111.sample000 = center_sample;

        // distibute faces left-right
        this.node100.sample111 = center_right_sample;
        this.node101.sample110 = center_right_sample;
        this.node110.sample101 = center_right_sample;
        this.node111.sample100 = center_right_sample;

        this.node000.sample011 = center_left_sample;
        this.node001.sample010 = center_left_sample;
        this.node010.sample001 = center_left_sample;
        this.node011.sample000 = center_left_sample;

        // distibute faces near-far
        this.node010.sample111 = center_far_sample;
        this.node011.sample110 = center_far_sample;
        this.node110.sample011 = center_far_sample;
        this.node111.sample010 = center_far_sample;

        this.node000.sample101 = center_near_sample;
        this.node001.sample100 = center_near_sample;
        this.node100.sample001 = center_near_sample;
        this.node101.sample000 = center_near_sample;

        // distibute faces bottom-top
        this.node001.sample111 = center_top_sample;
        this.node011.sample101 = center_top_sample;
        this.node101.sample011 = center_top_sample;
        this.node111.sample001 = center_top_sample;

        this.node000.sample110 = center_bottom_sample;
        this.node010.sample100 = center_bottom_sample;
        this.node100.sample010 = center_bottom_sample;
        this.node110.sample000 = center_bottom_sample;

        // distibute edges
        this.node100.sample000 = this.node000.sample100 = center_min_min_sample;
        this.node110.sample010 = this.node010.sample110 = center_max_min_sample;
        this.node111.sample011 = this.node011.sample111 = center_max_max_sample;
        this.node101.sample001 = this.node001.sample101 = center_min_max_sample;
        this.node010.sample000 = this.node000.sample010 = min_center_min_sample;
        this.node110.sample100 = this.node100.sample110 = max_center_min_sample;
        this.node111.sample101 = this.node101.sample111 = max_center_max_sample;
        this.node011.sample001 = this.node001.sample011 = min_center_max_sample;
        this.node001.sample000 = this.node000.sample001 = min_min_center_sample;
        this.node101.sample100 = this.node100.sample101 = max_min_center_sample;
        this.node111.sample110 = this.node110.sample111 = max_max_center_sample;
        this.node011.sample010 = this.node010.sample011 = min_max_center_sample;

        this.sample000 = null;
        this.sample001 = null;
        this.sample010 = null;
        this.sample011 = null;
        this.sample100 = null;
        this.sample101 = null;
        this.sample110 = null;
        this.sample111 = null;

        this.node000 = this.node000.sampleSignedDistanceField(field, merge_threshold);
        this.node001 = this.node001.sampleSignedDistanceField(field, merge_threshold);
        this.node010 = this.node010.sampleSignedDistanceField(field, merge_threshold);
        this.node011 = this.node011.sampleSignedDistanceField(field, merge_threshold);
        this.node100 = this.node100.sampleSignedDistanceField(field, merge_threshold);
        this.node101 = this.node101.sampleSignedDistanceField(field, merge_threshold);
        this.node110 = this.node110.sampleSignedDistanceField(field, merge_threshold);
        this.node111 = this.node111.sampleSignedDistanceField(field, merge_threshold);

        const children_kind = this.node000.kind | this.node100.kind | this.node110.kind | this.node010.kind
        |                     this.node001.kind | this.node101.kind | this.node111.kind | this.node011.kind;

        if (children_kind == SAMPLE_OCTREE_KIND_FULL) {
            return _SAMPLE_OCTREE_FULL;
        }

        if (children_kind == SAMPLE_OCTREE_KIND_EMPTY) {
            return _SAMPLE_OCTREE_EMPTY;
        }

        return this;
    }
}

export class SampleOctreeLeaf {
    constructor(min_x, min_y, min_z, max_x, max_y, max_z) {
        this.kind  = SAMPLE_OCTREE_KIND_LEAF;
        this.min_x = min_x;
        this.min_y = min_y;
        this.min_z = min_z;
        this.max_x = max_x;
        this.max_y = max_y;
        this.max_z = max_z;
        this.sample000 = 0;
        this.sample100 = 0;
        this.sample110 = 0;
        this.sample010 = 0;
        this.sample001 = 0;
        this.sample101 = 0;
        this.sample111 = 0;
        this.sample011 = 0;
        this.surface_vertex = null;
    }

    getNodeState() {
        let node_state = 0;
        if (this.sample000 < 0) node_state |= 0x01;
        if (this.sample001 < 0) node_state |= 0x02;
        if (this.sample010 < 0) node_state |= 0x04;
        if (this.sample011 < 0) node_state |= 0x08;
        if (this.sample100 < 0) node_state |= 0x10;
        if (this.sample101 < 0) node_state |= 0x20;
        if (this.sample110 < 0) node_state |= 0x40;
        if (this.sample111 < 0) node_state |= 0x80;
        return node_state;
    }

    sampleSignedDistanceField(field) {
        const min_x = this.min_x, max_x = this.max_x, center_x = (max_x + min_x) / 2;
        const min_y = this.min_y, max_y = this.max_y, center_y = (max_y + min_y) / 2;
        const min_z = this.min_z, max_z = this.max_z, center_z = (max_z + min_z) / 2;
        const node_state = this.getNodeState();

        if (node_state == ALL_SAMPLES_INSIDE) { return _SAMPLE_OCTREE_FULL; }
        if (node_state == ALL_SAMPLES_OUTSIDE) { return _SAMPLE_OCTREE_EMPTY; }

        const position = this.calculateSurfacePosition();
        const normal = field.calculateGradient(position.x, position.y, position.z).normalize();
        this.surface_vertex = new Vertex(position, normal);
        return this;
    }
        
    calculateSurfacePosition() {
        const edges = this.getTemporaryEdges();

        let x = 0, y = 0, z = 0, intersection_count = 0;
        for (let idx = 0; idx < edges.length; ++idx) {
            const edge = edges[idx];
            if ((edge.sample0 < 0) != (edge.sample1 < 0)) {
                const t = Math.abs(edge.sample0) / (Math.abs(edge.sample0) + Math.abs(edge.sample1));
                x += (edge.x1 - edge.x0) * t + edge.x0;
                y += (edge.y1 - edge.y0) * t + edge.y0;
                z += (edge.z1 - edge.z0) * t + edge.z0;
                intersection_count += 1;
            }
        }

        if (intersection_count == 0) {
            return null;
        }

        const factor = 1 / intersection_count;
        return new Vector3(x * factor, y * factor, z * factor);
    }
    
    getTemporaryEdges() {
        const min_x = this.min_x, min_y = this.min_y, min_z = this.min_z;
        const max_x = this.max_x, max_y = this.max_y, max_z = this.max_z;

        _temp_edges[0].set(min_x, min_y, min_z, this.sample000, max_x, min_y, min_z, this.sample100);
        _temp_edges[1].set(min_x, max_y, min_z, this.sample010, max_x, max_y, min_z, this.sample110);
        _temp_edges[2].set(min_x, min_y, max_z, this.sample001, max_x, min_y, max_z, this.sample101);
        _temp_edges[3].set(min_x, max_y, max_z, this.sample011, max_x, max_y, max_z, this.sample111);
        _temp_edges[4].set(min_x, min_y, min_z, this.sample000, min_x, max_y, min_z, this.sample010);
        _temp_edges[5].set(max_x, min_y, min_z, this.sample100, max_x, max_y, min_z, this.sample110);
        _temp_edges[6].set(min_x, min_y, max_z, this.sample001, min_x, max_y, max_z, this.sample011);
        _temp_edges[7].set(max_x, min_y, max_z, this.sample101, max_x, max_y, max_z, this.sample111);
        _temp_edges[8].set(min_x, min_y, min_z, this.sample000, min_x, min_y, max_z, this.sample001);
        _temp_edges[9].set(max_x, min_y, min_z, this.sample100, max_x, min_y, max_z, this.sample101);
        _temp_edges[10].set(max_x, max_y, min_z, this.sample110, max_x, max_y, max_z, this.sample111);
        _temp_edges[11].set(min_x, max_y, min_z, this.sample010, min_x, max_y, max_z, this.sample011);

        return _temp_edges;
    }
}

export class SampleOctreeEmpty {
    constructor() {
        this.kind = SAMPLE_OCTREE_KIND_EMPTY;
    }
}

export class SampleOctreeFull {
    constructor() {
        this.kind = SAMPLE_OCTREE_KIND_FULL;
    }
}

export class SampleOctreeMeshBuilder {
    constructor() {
        this.vertices = [];
        this.indices = [];
    }

    getMesh() {
        return new TriangleMesh(this.vertices, this.indices);
    }

    

    clear() {
        this.vertices = [];
        this.indices = [];
    }

    buildMesh(... nodes) {
        for (let idx = 0; idx < nodes.length; ++idx) {
            this.processCell(nodes[idx]);
        }

        return new TriangleMesh(this.vertices, this.indices);
    }

    addTriangle(node0, node1, node2, should_flip = false) {
        const vertex_count = this.vertices.length;
        
        if (should_flip) {
            this.indices.push(vertex_count + 2, vertex_count + 1, vertex_count + 0);
        }
        else {
            this.indices.push(vertex_count + 0, vertex_count + 1, vertex_count + 2);
        }

        this.vertices.push(
            node0.surface_vertex,
            node1.surface_vertex,
            node2.surface_vertex,
        );
    }

    addQuad(node0, node1, node2, node3, should_flip = false) {
        if (node0 == node1) return this.addTriangle(node0, node2, node3, should_flip);
        if (node0 == node2) throw new Error("Invalid quad winding!");
        if (node0 == node3) return this.addTriangle(node0, node1, node2, should_flip);
        if (node1 == node2) return this.addTriangle(node0, node1, node3, should_flip);
        if (node1 == node3) throw new Error("Invalid quad winding!");
        if (node2 == node3) return this.addTriangle(node0, node1, node2, should_flip);

        const vertex_count = this.vertices.length;
        
        if (should_flip) {
            this.indices.push(vertex_count + 3, vertex_count + 2, vertex_count + 0,
                              vertex_count + 2, vertex_count + 1, vertex_count + 0);
        }
        else {
            this.indices.push(vertex_count + 0, vertex_count + 1, vertex_count + 2,
                              vertex_count + 0, vertex_count + 2, vertex_count + 3);
        }

        this.vertices.push(
            node0.surface_vertex,
            node1.surface_vertex,
            node2.surface_vertex,
            node3.surface_vertex,
        );
    }

    processCell(node) {
        if (node.kind != SAMPLE_OCTREE_KIND_ROOT) {
            return;
        }

        this.processCell(node.node000);
        this.processCell(node.node001);
        this.processCell(node.node010);
        this.processCell(node.node011);
        this.processCell(node.node100);
        this.processCell(node.node101);
        this.processCell(node.node110);
        this.processCell(node.node111);

        this.processFaceX(node.node000, node.node100);
        this.processFaceX(node.node010, node.node110);
        this.processFaceX(node.node011, node.node111);
        this.processFaceX(node.node001, node.node101);

        this.processFaceY(node.node000, node.node010);
        this.processFaceY(node.node100, node.node110);
        this.processFaceY(node.node101, node.node111);
        this.processFaceY(node.node001, node.node011);

        this.processFaceZ(node.node000, node.node001);
        this.processFaceZ(node.node100, node.node101);
        this.processFaceZ(node.node110, node.node111);
        this.processFaceZ(node.node010, node.node011);

        this.processEdgeX(node.node000, node.node010, node.node011, node.node001);
        this.processEdgeX(node.node100, node.node110, node.node111, node.node101);

        this.processEdgeY(node.node000, node.node100, node.node101, node.node001);
        this.processEdgeY(node.node010, node.node110, node.node111, node.node011);

        this.processEdgeZ(node.node000, node.node100, node.node110, node.node010);
        this.processEdgeZ(node.node001, node.node101, node.node111, node.node011);
    }

    processFaceX(node0, node1) {
        const kind = node0.kind | node1.kind;
        if ((kind & SAMPLE_OCTREE_FULL_OR_EMPTY) || (kind & SAMPLE_OCTREE_KIND_ROOT) == 0) {
            return;
        }

        let children000 = node0;
        let children010 = node0;
        let children100 = node1;
        let children110 = node1;
        let children001 = node0;
        let children011 = node0;
        let children101 = node1;
        let children111 = node1;

        if (node0.kind == SAMPLE_OCTREE_KIND_ROOT) {
            children000 = node0.node100;
            children010 = node0.node110;
            children001 = node0.node101;
            children011 = node0.node111;
        }

        if (node1.kind == SAMPLE_OCTREE_KIND_ROOT) {
            children100 = node1.node000;
            children110 = node1.node010;
            children101 = node1.node001;
            children111 = node1.node011;
        }

        this.processFaceX(children000, children100);
        this.processFaceX(children010, children110);
        this.processFaceX(children011, children111);
        this.processFaceX(children001, children101);

        this.processEdgeY(children000, children100, children101, children001);
        this.processEdgeY(children010, children110, children111, children011);

        this.processEdgeZ(children000, children100, children110, children010);
        this.processEdgeZ(children001, children101, children111, children011);
    }

    processFaceY(node0, node1) {
        const kind = node0.kind | node1.kind;
        if ((kind & SAMPLE_OCTREE_FULL_OR_EMPTY) || (kind & SAMPLE_OCTREE_KIND_ROOT) == 0) {
            return;
        }

        let children000 = node0;
        let children100 = node0;
        let children010 = node1;
        let children110 = node1;
        let children001 = node0;
        let children101 = node0;
        let children011 = node1;
        let children111 = node1;

        if (node0.kind == SAMPLE_OCTREE_KIND_ROOT) {
            children000 = node0.node010;
            children100 = node0.node110;
            children001 = node0.node011;
            children101 = node0.node111;
        }

        if (node1.kind == SAMPLE_OCTREE_KIND_ROOT) {
            children010 = node1.node000;
            children110 = node1.node100;
            children011 = node1.node001;
            children111 = node1.node101;
        }

        this.processFaceY(children000, children010);
        this.processFaceY(children100, children110);
        this.processFaceY(children101, children111);
        this.processFaceY(children001, children011);

        this.processEdgeX(children000, children010, children011, children001);
        this.processEdgeX(children100, children110, children111, children101);

        this.processEdgeZ(children000, children100, children110, children010);
        this.processEdgeZ(children001, children101, children111, children011);
    }

    processFaceZ(node0, node1) {
        const kind = node0.kind | node1.kind;
        if ((kind & SAMPLE_OCTREE_FULL_OR_EMPTY) || (kind & SAMPLE_OCTREE_KIND_ROOT) == 0) {
            return;
        }

        let children000 = node0;
        let children100 = node0;
        let children001 = node1;
        let children101 = node1;
        let children010 = node0;
        let children110 = node0;
        let children011 = node1;
        let children111 = node1;

        if (node0.kind == SAMPLE_OCTREE_KIND_ROOT) {
            children000 = node0.node001;
            children100 = node0.node101;
            children010 = node0.node011;
            children110 = node0.node111;
        }

        if (node1.kind == SAMPLE_OCTREE_KIND_ROOT) {
            children001 = node1.node000;
            children101 = node1.node100;
            children011 = node1.node010;
            children111 = node1.node110;
        }

        this.processFaceZ(children000, children001);
        this.processFaceZ(children100, children101);
        this.processFaceZ(children110, children111);
        this.processFaceZ(children010, children011);

        this.processEdgeX(children000, children010, children011, children001);
        this.processEdgeX(children100, children110, children111, children101);

        this.processEdgeY(children000, children100, children101, children001);
        this.processEdgeY(children010, children110, children111, children011);
    }

    processEdgeX(node0, node1, node2, node3) {
        const kind = node0.kind | node1.kind | node2.kind | node3.kind;

        if (kind & SAMPLE_OCTREE_FULL_OR_EMPTY) {
            return;
        }

        if ((kind & SAMPLE_OCTREE_KIND_ROOT) == 0) {
            const test = getLowestSubdivision(node0, node1, node2, node3);

            let should_flip = false;
            switch (test) {
            case node0: should_flip = node0.sample111 < node0.sample011; break;
            case node1: should_flip = node1.sample101 < node1.sample001; break;
            case node2: should_flip = node2.sample100 < node2.sample000; break;
            case node3: should_flip = node3.sample110 < node3.sample010; break;
            }

            this.addQuad(node0, node1, node2, node3, should_flip);
        }
        else {
            const children000 = node0.kind == SAMPLE_OCTREE_KIND_ROOT ? node0.node011 : node0;
            const children100 = node0.kind == SAMPLE_OCTREE_KIND_ROOT ? node0.node111 : node0;
            const children010 = node1.kind == SAMPLE_OCTREE_KIND_ROOT ? node1.node001 : node1;
            const children110 = node1.kind == SAMPLE_OCTREE_KIND_ROOT ? node1.node101 : node1;
            const children001 = node3.kind == SAMPLE_OCTREE_KIND_ROOT ? node3.node010 : node3;
            const children101 = node3.kind == SAMPLE_OCTREE_KIND_ROOT ? node3.node110 : node3;
            const children011 = node2.kind == SAMPLE_OCTREE_KIND_ROOT ? node2.node000 : node2;
            const children111 = node2.kind == SAMPLE_OCTREE_KIND_ROOT ? node2.node100 : node2;

            this.processEdgeX(children000, children010, children011, children001);
            this.processEdgeX(children100, children110, children111, children101);
        }
    }

    processEdgeY(node0, node1, node2, node3) {
        const kind = node0.kind | node1.kind | node2.kind | node3.kind;

        if (kind & SAMPLE_OCTREE_FULL_OR_EMPTY) {
            return;
        }

        if ((kind & SAMPLE_OCTREE_KIND_ROOT) == 0) {
            const test = getLowestSubdivision(node0, node1, node2, node3);

            let should_flip = false;
            switch (test) {
            case node0: should_flip = node0.sample101 < node0.sample111; break;
            case node1: should_flip = node1.sample001 < node1.sample011; break;
            case node2: should_flip = node2.sample000 < node2.sample010; break;
            case node3: should_flip = node3.sample100 < node3.sample110; break;
            }

            this.addQuad(node0, node1, node2, node3, should_flip);
        }
        else {
            const children000 = node0.kind == SAMPLE_OCTREE_KIND_ROOT ? node0.node101 : node0;
            const children010 = node0.kind == SAMPLE_OCTREE_KIND_ROOT ? node0.node111 : node0;
            const children100 = node1.kind == SAMPLE_OCTREE_KIND_ROOT ? node1.node001 : node1;
            const children110 = node1.kind == SAMPLE_OCTREE_KIND_ROOT ? node1.node011 : node1;
            const children001 = node3.kind == SAMPLE_OCTREE_KIND_ROOT ? node3.node100 : node3;
            const children011 = node3.kind == SAMPLE_OCTREE_KIND_ROOT ? node3.node110 : node3;
            const children101 = node2.kind == SAMPLE_OCTREE_KIND_ROOT ? node2.node000 : node2;
            const children111 = node2.kind == SAMPLE_OCTREE_KIND_ROOT ? node2.node010 : node2;

            this.processEdgeY(children000, children100, children101, children001);
            this.processEdgeY(children010, children110, children111, children011);
        }
    }

    processEdgeZ(node0, node1, node2, node3) {
        const kind = node0.kind | node1.kind | node2.kind | node3.kind;

        if (kind & SAMPLE_OCTREE_FULL_OR_EMPTY) {
            return;
        }

        if ((kind & SAMPLE_OCTREE_KIND_ROOT) == 0) {
            const test = getLowestSubdivision(node0, node1, node2, node3);

            let should_flip = false;
            switch (test) {
            case node0: should_flip = node0.sample111 < node0.sample110; break;
            case node1: should_flip = node1.sample011 < node1.sample010; break;
            case node2: should_flip = node2.sample001 < node2.sample000; break;
            case node3: should_flip = node3.sample101 < node3.sample100; break;
            }

            this.addQuad(node0, node1, node2, node3, should_flip);
        }
        else {
            const children000 = node0.kind == SAMPLE_OCTREE_KIND_ROOT ? node0.node110 : node0;
            const children001 = node0.kind == SAMPLE_OCTREE_KIND_ROOT ? node0.node111 : node0;
            const children100 = node1.kind == SAMPLE_OCTREE_KIND_ROOT ? node1.node010 : node1;
            const children101 = node1.kind == SAMPLE_OCTREE_KIND_ROOT ? node1.node011 : node1;
            const children010 = node3.kind == SAMPLE_OCTREE_KIND_ROOT ? node3.node100 : node3;
            const children011 = node3.kind == SAMPLE_OCTREE_KIND_ROOT ? node3.node101 : node3;
            const children110 = node2.kind == SAMPLE_OCTREE_KIND_ROOT ? node2.node000 : node2;
            const children111 = node2.kind == SAMPLE_OCTREE_KIND_ROOT ? node2.node001 : node2;

            this.processEdgeZ(children000, children100, children110, children010);
            this.processEdgeZ(children001, children101, children111, children011);
        }
    }
}

export function createSampleOctreeFromField(field, bounding_box, subdivision_count, merge_threshold = null) {
    return new SampleOctreeRoot(
        bounding_box.min_x,
        bounding_box.min_y,
        bounding_box.min_z,
        bounding_box.max_x,
        bounding_box.max_y,
        bounding_box.max_z,
        subdivision_count
    ).sampleSignedDistanceField(field, merge_threshold);
}

export function triangulateSignedDistanceField(field, bounding_box = field.bounding_box, subdivision_count = 8, merge_threshold = 0.0005) {
    const tree = new SampleOctreeRoot(
        bounding_box.min_x,
        bounding_box.min_y,
        bounding_box.min_z,
        bounding_box.max_x,
        bounding_box.max_y,
        bounding_box.max_z,
        subdivision_count
    ).sampleSignedDistanceField(field, merge_threshold);
    const _builder = new SampleOctreeMeshBuilder();
    _builder.processCell(tree);
    return _builder.getMesh();
}

export class BufferWriter {
    constructor(buffer) {
        this.view = new DataView(buffer);
        this.offset = 0;
    }

    setOffset(offset) { this.offset = offset; }
    writeInt8(value) { this.view.setInt8(this.offset, value); this.offset += 1; }
    writeUint8(value) { this.view.setUint8(this.offset, value); this.offset += 1; }
    writeInt16(value, littleEndian) { this.view.setInt16(this.offset, value, littleEndian); this.offset += 2; }
    writeUint16(value, littleEndian) { this.view.setUint16(this.offset, value, littleEndian); this.offset += 2; }
    writeInt32(value, littleEndian) { this.view.setInt32(this.offset, value, littleEndian); this.offset += 4; }
    writeUint32(value, littleEndian) { this.view.setUint32(this.offset, value, littleEndian); this.offset += 4; }
    writeFloat32(value, littleEndian) { this.view.setFloat32(this.offset, value, littleEndian); this.offset += 4; }
    writeFloat64(value, littleEndian) { this.view.setFloat64(this.offset, value, littleEndian); this.offset += 8; }
    writeBuffer(buffer) {
        const uint8_view = new Uint8Array(this.view.buffer);
        uint8_view.set(new Uint8Array(buffer), this.offset);
        this.offset += buffer.byteLength;
    }
}

export function convertToSTL(mesh) {
    if (mesh instanceof TriangleMesh == false) {
        throw new Error("input mesh must be a 'TriangleMesh'!");
    }

    const triangle_count = mesh.indices.length / 3;
    const buffer = new ArrayBuffer(STL_TRIANGLE_OFFSET + STL_BYTES_PER_TRIANGLE * triangle_count);
    const writer = new BufferWriter(buffer);

    writer.setOffset(STL_HEADER_SIZE);
    writer.writeUint32(triangle_count, true);

    for (let idx = 0; idx < mesh.indices.length; idx += 3) {
        const position0 = mesh.vertices[mesh.indices[idx + 0]].position;
        const position1 = mesh.vertices[mesh.indices[idx + 1]].position;
        const position2 = mesh.vertices[mesh.indices[idx + 2]].position;

        temp_plane.fromVectors(position0, position1, position2);

        writer.writeFloat32(temp_plane.normal_x, true);
        writer.writeFloat32(temp_plane.normal_y, true);
        writer.writeFloat32(temp_plane.normal_z, true);

        writer.writeFloat32(position0.x, true);
        writer.writeFloat32(position0.y, true);
        writer.writeFloat32(position0.z, true);

        writer.writeFloat32(position1.x, true);
        writer.writeFloat32(position1.y, true);
        writer.writeFloat32(position1.z, true);

        writer.writeFloat32(position2.x, true);
        writer.writeFloat32(position2.y, true);
        writer.writeFloat32(position2.z, true);

        writer.writeUint16(0, true);
    }

    return new Uint8Array(buffer);
}

export function convertToOBJ(mesh) {
    if (mesh instanceof TriangleMesh == false) {
        throw new Error("input mesh must be a 'TriangleMesh'!");
    }

    let output = "";
    for (let idx = 0; idx < mesh.vertices.length; ++idx) {
        const position = mesh.vertices[idx].position;
        output += `v ${position.x}, ${position.y}, ${position.z}\n`;
    }

    for (let idx = 0; idx < mesh.indices.length; idx += 3) {
        const index0 = mesh.indices[idx + 0] + 1;
        const index1 = mesh.indices[idx + 1] + 1;
        const index2 = mesh.indices[idx + 2] + 1;
        output += `f ${index0}, ${index1}, ${index2}\n`;
    }

    return output;
}

export class SignedDistanceField3 {
    calculateSignedDistance(x, y, z) {
        throw new Error("calculateSignedDistance is unimplemented!");
    }

    calculateGradient(x, y, z, out_vector = new Vector3()) {
        const EPSILON = 0.00001;
        const sample0 = this.calculateSignedDistance(x + EPSILON, y - EPSILON, z - EPSILON);
        const sample1 = this.calculateSignedDistance(x - EPSILON, y - EPSILON, z + EPSILON);
        const sample2 = this.calculateSignedDistance(x - EPSILON, y + EPSILON, z - EPSILON);
        const sample3 = this.calculateSignedDistance(x + EPSILON, y + EPSILON, z + EPSILON);
        return out_vector.set(
            + sample0 - sample1 - sample2 + sample3,
            - sample0 - sample1 + sample2 + sample3,
            - sample0 + sample1 - sample2 + sample3,
        ).normalize();
    }

    raycast(x0, y0, z0, dx, dy, dz, distance_max = Infinity) {
        let total_distance = 0;
        for (let idx = 0; idx < INTERSECTION_DEPTH_MAX && total_distance < distance_max; ++idx) {
            const radius = this.calculateSignedDistance(
                x0 + total_distance * dx,
                y0 + total_distance * dy,
                z0 + total_distance * dz);

            if (radius < EPSILON) {
                return total_distance;
            }
            
            total_distance += radius;
        }

        return null;
    }

    raycastIntersection(x0, y0, z0, x1, y1, z1, out = new Vector3()) {
        const sample0 = this.calculateSignedDistance(x0, y0, z0);
        const sample1 = this.calculateSignedDistance(x1, y1, z1);

        if ((sample0 < 0) == (sample1 < 0)) {
            return null;
        }

        if (sample0 < 0) {
            [ x0, y0, z0, x1, y1, z1 ] = [ x1, y1, z1, x0, y0, z0 ];
        }

        const dx = x1 - x0;
        const dy = y1 - y0;
        const dz = z1 - z0;

        const length = Math.sqrt(dx*dx + dy*dy + dz*dz)
        
        if (length == 0) {
            return null;
        }

        const inv = 1 / length;
        const normal_x = dx * inv;
        const normal_y = dy * inv;
        const normal_z = dz * inv;

        const distance = this.raycast(x0, y0, z0, normal_x, normal_y, normal_z, length);

        if (distance != null) {
            return out.set(x0 + distance * normal_x,
                           y0 + distance * normal_y,
                           z0 + distance * normal_z);
        }

        return null;
    }

    _raycastIntersection(x0, y0, z0, x1, y1, z1, out = new Vector3()) {
        const dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
        const inv = 1 / Math.sqrt(dx*dx + dy*dy + dz*dz);
        const nx = dx * inv, ny = dy * inv, nz = dz * inv;

        for (let idx = 0; idx < INTERSECTION_DEPTH_MAX; ++idx) {
            const sphere_distance = this.calculateSignedDistance(x0, y0, z0);
            const positive_sphere_distance = Math.abs(sphere_distance);

            x0 += nx * positive_sphere_distance;
            y0 += ny * positive_sphere_distance;
            z0 += nz * positive_sphere_distance;
        }

        return null;
    }

    transform(matrix = new Matrix4()) { return new TransformSignedDistanceField3(this, matrix); }
    rotateX(angle = 0) { return this.transform(_m4_temp.identity().rotateX(angle)); }
    rotateY(angle = 0) { return this.transform(_m4_temp.identity().rotateY(angle)); }
    rotateZ(angle = 0) { return this.transform(_m4_temp.identity().rotateZ(angle)); }
    twist(angle = 0) { return new TwistSignedDistanceField3(this, angle); }
    shell(thickness = 1) { return new ShellSignedDistanceField3(this, thickness); }
    offset(distance = 1) { return new OffsetSignedDistanceField3(this, distance); }
    union(... fields) { return new UnionSignedDistanceField3([ this, ... fields ]); }
    difference(... fields) { return new DifferenceSignedDistanceField3([ this, ... fields ]); }
    intersection(... fields) { return new IntersectionSignedDistanceField3([ this, ... fields ]); }

    translate(x = 0, y = 0, z = 0) {
        return this.transform(_m4_temp.identity().translate(x, y, z));
    }

    rotate(x = 0, y = 0, z = 1, angle = 0) {
        return this.transform(_m4_temp.identity().rotate(x, y, z, angle));
    }
}

export class BoxSignedDistanceField3 extends SignedDistanceField3 {
    constructor(size_x, size_y, size_z, corner_radius = 0) {
        super();
        this.corner_radius = corner_radius;
        this.half_size_x = 0.5 * size_x;
        this.half_size_y = 0.5 * size_y;
        this.half_size_z = 0.5 * size_z;
        this.bounding_box = new BoundingBox3(-this.half_size_x, -this.half_size_y, -this.half_size_z,
                                              this.half_size_x,  this.half_size_y,  this.half_size_z);
    }

    calculateSignedDistance(x, y, z) {
        const pos_x = Math.abs(x) - this.half_size_x + this.corner_radius;
        const pos_y = Math.abs(y) - this.half_size_y + this.corner_radius;
        const pos_z = Math.abs(z) - this.half_size_z + this.corner_radius;
        const dx = Math.max(pos_x, 0);
        const dy = Math.max(pos_y, 0);
        const dz = Math.max(pos_z, 0);
        const side = Math.min(Math.max(pos_x, pos_y, pos_z), 0.0);
        return Math.sqrt(dx*dx + dy*dy + dz*dz) + side - this.corner_radius;
    }
}

export class BoxFrameSignedDistanceField3 extends SignedDistanceField3 {
    constructor(size_x, size_y, size_z, frame_size = 5) {
        super();
        this.frame_size = frame_size;
        this.half_size_x = 0.5 * size_x;
        this.half_size_y = 0.5 * size_y;
        this.half_size_z = 0.5 * size_z;
        this.bounding_box = new BoundingBox3(-this.half_size_x, -this.half_size_y, -this.half_size_z,
                                              this.half_size_x,  this.half_size_y,  this.half_size_z);
    }

    calculateSignedDistance(x, y, z) {
        const p0 = Math.abs(x) - this.half_size_x;
        const p1 = Math.abs(y) - this.half_size_y;
        const p2 = Math.abs(z) - this.half_size_z;
        const q0 = Math.abs(p0 + this.frame_size) - this.frame_size;
        const q1 = Math.abs(p1 + this.frame_size) - this.frame_size;
        const q2 = Math.abs(p2 + this.frame_size) - this.frame_size;

        return Math.min(
            distanceFromOrigin(Math.max(q0, 0), Math.max(q1, 0), Math.max(p2, 0)) + Math.min(Math.max(q0, q1, p2), 0.0),
            distanceFromOrigin(Math.max(p0, 0), Math.max(q1, 0), Math.max(q2, 0)) + Math.min(Math.max(p0, q1, q2), 0.0),
            distanceFromOrigin(Math.max(q0, 0), Math.max(p1, 0), Math.max(q2, 0)) + Math.min(Math.max(q0, p1, q2), 0.0)
        );
    }
}

export class SphereSignedDistanceField3 extends SignedDistanceField3 {
    constructor(radius = 5) {
        super();
        this.radius = radius;
        this.bounding_box = new BoundingBox3(-radius, -radius, -radius,
                                              radius,  radius,  radius);
    }

    calculateSignedDistance(x, y, z) {
        return Math.sqrt(x*x + y*y + z*z) - this.radius;
    }
}

export class TorusSignedDistanceField3 extends SignedDistanceField3 {
    constructor(outer_radius = 10, inner_radius = 5) {
        super();
        const farthest_radius = outer_radius + inner_radius;
        this.inner_radius = inner_radius;
        this.outer_radius = outer_radius;
        this.bounding_box = new BoundingBox3(-farthest_radius, -farthest_radius, -inner_radius,
                                              farthest_radius,  farthest_radius,  inner_radius);
    }

    calculateSignedDistance(x, y, z) {
        const distance_xy = Math.sqrt(x*x + y*y) - this.outer_radius;
        return Math.sqrt(distance_xy * distance_xy + z*z) - this.inner_radius;
    }
}

export class CylinderSignedDistanceField3 extends SignedDistanceField3 {
    constructor(radius = 10, height = 20) {
        super();
        this.radius = radius;
        this.half_height = 0.5 * height;
        this.bounding_box = new BoundingBox3(-this.radius, -this.radius, -this.half_height,
                                              this.radius,  this.radius,  this.half_height);
    }

    calculateSignedDistance(x, y, z) {
        const d0 = Math.abs(Math.sqrt(x*x + y*y)) - this.radius;
        const d1 = Math.abs(z) - this.half_height;
        const a0 = Math.max(d0,0.0);
        const a1 = Math.max(d1,0.0);
        return Math.min(Math.max(d0,d1),0.0) + Math.sqrt(a0*a0 + a1*a1);
    }
}

export class PlaneSignedDistanceField3 extends SignedDistanceField3 {
    constructor(plane = new Plane(0, 0, 1, 0)) {
        super();
        this.plane = plane;
        this.bounding_box = new BoundingBox3(Infinity,  Infinity,  Infinity,
                                            -Infinity, -Infinity, -Infinity);
    }

    calculateSignedDistance(x, y, z) {
        return this.plane.distanceToPoint(x, y, z);
    }
}

export class TransformSignedDistanceField3 extends SignedDistanceField3 {
    constructor(field, matrix = new Matrix4()) {
        super();
        this.matrix = matrix.clone();
        this.inverse_matrix = this.matrix.clone().inverse();
        this.bounding_box = field.bounding_box.clone().transformMatrix4(this.matrix);
        this.field = field;
    }

    calculateSignedDistance(x, y, z) {
        _v3_temp.set(x, y, z).transformMatrix4(this.inverse_matrix);
        return this.field.calculateSignedDistance(_v3_temp.x, _v3_temp.y, _v3_temp.z);
    }
}

export class TwistSignedDistanceField3 extends SignedDistanceField3 {
    constructor(field, angle) {
        super();
        this.field = field;
        this.k = angle / (field.bounding_box.max_z - field.bounding_box.min_z);
        const min = Math.min(field.bounding_box.min_x, field.bounding_box.min_y);
        const max = Math.max(field.bounding_box.max_x, field.bounding_box.max_y);
        this.bounding_box = new BoundingBox3(min, min, field.bounding_box.min_z,
                                             max, max, field.bounding_box.max_z);
    }

    calculateSignedDistance(x, y, z) {
        const c = Math.cos(z * this.k);
        const s = Math.sin(z * this.k);
        const rotated_x = x * c - y * s;
        const rotated_y = x * s + y * c;
        return this.field.calculateSignedDistance(rotated_x, rotated_y, z);
    }
}

export class ShellSignedDistanceField3 extends SignedDistanceField3 {
    constructor(field, thickness) {
        super();
        this.field = field;
        this.half_thickness = 0.5 * thickness;
        this.bounding_box = field.bounding_box.clone().grow(this.half_thickness);
    }

    calculateSignedDistance(x, y, z) {
        return Math.abs(this.field.calculateSignedDistance(x, y, z)) - this.half_thickness;
    }
}

export class OffsetSignedDistanceField3 extends SignedDistanceField3 {
    constructor(field, distance) {
        super();
        this.field = field;
        this.distance = distance;
        this.bounding_box = field.bounding_box.clone().grow(this.distance);
    }

    calculateSignedDistance(x, y, z) {
        return this.field.calculateSignedDistance(x, y, z) - this.distance;
    }
}

export class UnionSignedDistanceField3 extends SignedDistanceField3 {
    constructor(fields = []) {
        super();
        this.fields = fields;
        this.bounding_box = this.fields[0].bounding_box.clone();
        for (let idx = 1; idx < this.fields.length; ++idx) {
            this.bounding_box.union(this.fields[idx].bounding_box);
        }
    }

    calculateSignedDistance(x, y, z) {
        let distance = this.fields[0].calculateSignedDistance(x, y, z);
        for (let idx = 1; idx < this.fields.length; ++idx) {
            distance = Math.min(distance, this.fields[idx].calculateSignedDistance(x, y, z));
        }
        return distance;
    }
}

export class DifferenceSignedDistanceField3 extends SignedDistanceField3 {
    constructor(fields = []) {
        super();
        this.fields = fields;
        this.bounding_box = this.fields[0].bounding_box.clone();
    }

    calculateSignedDistance(x, y, z) {
        let distance = this.fields[0].calculateSignedDistance(x, y, z);
        for (let idx = 1; idx < this.fields.length; ++idx) {
            distance = Math.max(distance, -this.fields[idx].calculateSignedDistance(x, y, z));
        }
        return distance;
    }
}

export class IntersectionSignedDistanceField3 extends SignedDistanceField3 {
    constructor(fields = []) {
        super();
        this.fields = fields;
        this.bounding_box = this.fields[0].bounding_box.clone();
        for (let idx = 1; idx < this.fields.length; ++idx) {
            this.bounding_box.intersect(this.fields[idx].bounding_box);
        }
    }

    calculateSignedDistance(x, y, z) {
        let distance = this.fields[0].calculateSignedDistance(x, y, z);
        for (let idx = 1; idx < this.fields.length; ++idx) {
            distance = Math.max(distance, this.fields[idx].calculateSignedDistance(x, y, z));
        }
        return distance;
    }
}


export class UnionSmoothSignedDistanceField3 extends SignedDistanceField3 {
    constructor(smoothness, ... fields) {
        super();
        this.smoothness = smoothness;
        this.fields = fields;
        this.bounding_box = this.fields[0].bounding_box.clone();
        for (let idx = 1; idx < this.fields.length; ++idx) {
            this.bounding_box.union(this.fields[idx].bounding_box);
        }
    }

    calculateSignedDistance(x, y, z) {
        let distance0 = this.fields[0].calculateSignedDistance(x, y, z);

        for (let idx = 1; idx < this.fields.length; ++idx) {
            let distance1 = this.fields[idx].calculateSignedDistance(x, y, z);
            let h = clamp(0.5 + 0.5 * (distance1 - distance0) / this.smoothness, 0, 1);
            let m = distance1 + (distance0 - distance1) * h;
            distance0 = m - this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldDifferenceSmooth extends SignedDistanceField {
    constructor(smoothness, ... fields) {
        super();
        this.smoothness = smoothness;
        this.fields = fields;
        this.bounding_box = this.fields[0].bounding_box.clone();
    }

    calculateSignedDistance(x, y, z) {
        let distance0 = this.fields[0].calculateSignedDistance(x, y, z);

        for (let idx = 1; idx < this.fields.length; ++idx) {
            let distance1 = this.fields[idx].calculateSignedDistance(x, y, z);
            let h = clamp(0.5 - 0.5 * (distance1 + distance0) / this.smoothness, 0, 1);
            let m = distance0 + (-distance1 - distance0) * h;
            distance0 = m + this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldIntersectionSmooth extends SignedDistanceField {
    constructor(smoothness, ... fields) {
        super();
        this.smoothness = smoothness;
        this.fields = fields;
        this.bounding_box = this.fields[0].bounding_box.clone();
        for (let idx = 1; idx < this.fields.length; ++idx) {
            this.bounding_box.intersect(this.fields[idx].bounding_box);
        }
    }

    calculateSignedDistance(x, y, z) {
        let distance0 = this.fields[0].calculateSignedDistance(x, y, z);

        for (let idx = 1; idx < this.fields.length; ++idx) {
            let distance1 = this.fields[idx].calculateSignedDistance(x, y, z);
            let h = clamp(0.5 - 0.5 * (distance1 - distance0) / this.smoothness, 0, 1);
            let m = distance1 + (distance0 - distance1) * h;
            distance0 = m + this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export function box(size_x = 20, size_y = 20, size_z = 20, corner_radius = 0) { return new BoxSignedDistanceField3(size_x, size_y, size_z, corner_radius); }
export function boxFrame(size_x = 20, size_y = 20, size_z = 20, frame_size = 5) { return new BoxFrameSignedDistanceField3(size_x, size_y, size_z, frame_size); }
export function sphere(radius = 5) { return new SphereSignedDistanceField3(radius); }
export function torus(radius_outer = 20, radius_inner = 5) { return new TorusSignedDistanceField3(radius_outer, radius_inner); }
export function cylinder(radius = 10, height = Infinity) { return new CylinderSignedDistanceField3(radius, height); }
export function plane(plane_ = new Plane()) { return new PlaneSignedDistanceField3(plane_); }
export function transform(field, matrix = new Matrix4()) { return new TransformSignedDistanceField3(field, matrix); }
export function twist(field, angle) { return new TwistSignedDistanceField3(field, angle); }
export function shell(field, thickness) { return new ShellSignedDistanceField3(field, thickness); }
export function offset(field, distance) { return new OffsetSignedDistanceField3(field, distance); }
export function union(... fields) { return new UnionSignedDistanceField3(fields); }
export function difference(... fields) { return new DifferenceSignedDistanceField3(fields); }
export function intersection(... fields) { return new IntersectionSignedDistanceField3(fields); }

const _SAMPLE_OCTREE_FULL = new SampleOctreeFull();
const _SAMPLE_OCTREE_EMPTY = new SampleOctreeEmpty();
const temp_plane = new Plane();
const _m4_temp = new Matrix4();
const _v3_temp = new Vector3(0, 0, 0);

const _temp_edges = [
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
    new SampleEdge(0, 0, 0, 0, 0, 0, 0, 0),
];
