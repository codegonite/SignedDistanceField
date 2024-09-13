// https://shaderfun.com/2018/07/23/signed-distance-fields-part-8-gradients-bevels-and-noise/
// https://www.boristhebrave.com/2018/04/15/dual-contouring-tutorial/
// https://github.com/ssloy/least-squares-course
// https://www.mattkeeter.com/projects/contours/
// https://iquilezles.org/articles/distfunctions/

/**
 * @template T
 * @typedef {[number, number, number]} Number3
 * @typedef {{x: number, y: number, z: number}} VectorObject3
 * @typedef {[number, number, number, number]} Number4
 * @typedef {{x: number, y: number, z: number, w: number}} VectorObject4
 * @typedef {[number, number, number, number,
 *            number, number, number, number,
 *            number, number, number, number,
 *            number, number, number, number]} Number16
 * @typedef {{m00: number, m10: number, m20: number, m30: number,
 *            m01: number, m11: number, m21: number, m31: number,
 *            m02: number, m12: number, m22: number, m32: number,
 *            m03: number, m13: number, m23: number, m33: number}} MatrixObject4
 */

export const EPSILON         = 1e-6;
export const INV_TWO_EPSILON = 1 / (2 * EPSILON);

export const SIGNED_DISTANCE_FIELD_KIND_UNKNOWN             = 0;
export const SIGNED_DISTANCE_FIELD_KIND_SPHERE              = 1;
export const SIGNED_DISTANCE_FIELD_KIND_CYLINDER            = 2;
export const SIGNED_DISTANCE_FIELD_KIND_CUBE                = 3;
export const SIGNED_DISTANCE_FIELD_KIND_BOX              = 4;
export const SIGNED_DISTANCE_FIELD_KIND_UNION               = 5;
export const SIGNED_DISTANCE_FIELD_KIND_DIFFERENCE          = 6;
export const SIGNED_DISTANCE_FIELD_KIND_INTERSECTION        = 7;
export const SIGNED_DISTANCE_FIELD_KIND_UNION_SMOOTH        = 8;
export const SIGNED_DISTANCE_FIELD_KIND_DIFFERENCE_SMOOTH   = 9;
export const SIGNED_DISTANCE_FIELD_KIND_INTERSECTION_SMOOTH = 10;
export const SIGNED_DISTANCE_FIELD_KIND_TRANSLATE           = 11;
export const SIGNED_DISTANCE_FIELD_KIND_TRANSFORM           = 12;
export const SIGNED_DISTANCE_FIELD_KIND_OFFSET              = 13;

/**
 * @param {number} value 
 * @param {number} min 
 * @param {number} max 
 * @returns {number}
 */
function clamp(value, min, max) {
    if (value > max) return max;
    if (value < min) return min;
    return value;
}

/**
 * @param {number} a
 * @param {number} b
 * @returns {number}
 */
function greatest_common_factor(a, b) {
    if (a < b) {
        return greatest_common_factor(b, a);
    }

    if (Math.abs(b) < 0.001) {
        return a;
    }

    return greatest_common_factor(b, a - Math.floor(a / b) * b);
}

export class Vector3 {
    /**
     * @param {VectorObject3} v
     * @returns {Vector3}
     */
    static from(v) {
        return new Vector3(v.x, v.y, v.z);
    }

    /**
     * @param {Number3} v
     * @returns {Vector3}
     */
    static fromArray(v) {
        return new Vector3(v[0], v[1], v[2]);
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @constructor
     */
    constructor(x=0, y=0, z=0) {
        /** @type {number} */
        this.x = x;

        /** @type {number} */
        this.y = y;

        /** @type {number} */
        this.z = z;
    }

    *[Symbol.iterator]() {
        yield this.x;
        yield this.y;
        yield this.z;
    }

    /** @returns {Vector3} */
    clone() { return new Vector3(this.x, this.y, this.z); }

    /** @returns {VectorObject3} */
    toObject() { return { x: this.x, y: this.y, z: this.z }; }

    /** @returns {Number3} */
    toArray() { return [ this.x, this.y, this.z ]; }

    /**
     * @param {VectorObject3} other
     * @returns {number}
     */
    distance(other) {
        const x = other.x - this.x;
        const y = other.y - this.y;
        const z = other.z - this.z;
        return Math.sqrt(x*x + y*y + z*z);
    }

    /**
     * @returns {number}
     */
    magnitudeSquared() { return this.x*this.x + this.y*this.y + this.z*this.z; }

    /**
     * @returns {number}
     */
    magnitude() { return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z); }

    /**
     * @param {VectorObject3} other
     * @returns {number}
     */
    dot(other) { return this.x*other.x + this.y*other.y + this.z*other.z; }

    /**
     * @param {VectorObject3} other
     * @param {number} t
     * @returns {Vector3}
     */
    lerp(other, t = 0) {
        const x0 = this.x,  y0 = this.y,  z0 = this.z;
        const x1 = other.x, y1 = other.y, z1 = other.z;

        return new Vector3((x1 - x0) * t + x0,
                           (y1 - y0) * t + y0,
                           (z1 - z0) * t + z0);
    }

    /**
     * @returns {Vector3}
     */
    abs() {
        this.x = Math.abs(this.x);
        this.y = Math.abs(this.y);
        this.z = Math.abs(this.z);
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    min(other) {
        this.x = Math.min(this.x, other.x);
        this.y = Math.min(this.y, other.y);
        this.z = Math.min(this.z, other.z);
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    max(other) {
        this.x = Math.max(this.x, other.x);
        this.y = Math.max(this.y, other.y);
        this.z = Math.max(this.z, other.z);
        return this;
    }

    /**
     * @param {VectorObject3} v0
     * @param {VectorObject3} v1
     * @returns {Vector3}
     */
    clamp(v0, v1) {
        this.x = clamp(this.x, v0.x, v1.x);
        this.y = clamp(this.y, v0.y, v1.y);
        this.z = clamp(this.z, v0.z, v1.z);
        return this;
    }

    /**
     * @returns {Vector3}
     */
    neg() {
        this.x = -this.x;
        this.y = -this.y;
        this.z = -this.z;
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    add(other) {
        this.x += other.x;
        this.y += other.y;
        this.z += other.z;
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    sub(other) {
        this.x -= other.x;
        this.y -= other.y;
        this.z -= other.z;
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    mul(other) {
        this.x *= other.x;
        this.y *= other.y;
        this.z *= other.z;
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    div(other) {
        this.x /= other.x;
        this.y /= other.y;
        this.z /= other.z;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector3}
     */
    addScalar(scalar) {
        this.x += scalar;
        this.y += scalar;
        this.z += scalar;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector3}
     */
    subScalar(scalar) {
        this.x -= scalar;
        this.y -= scalar;
        this.z -= scalar;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector3}
     */
    mulScalar(scalar) {
        this.x *= scalar;
        this.y *= scalar;
        this.z *= scalar;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector3}
     */
    divScalar(scalar) {
        this.x /= scalar;
        this.y /= scalar;
        this.z /= scalar;
        return this;
    }

    /**
     * @param {number} length
     * @returns {Vector3}
     */
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

    /**
     * @returns {Vector3}
     */
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

    /**
     * @param {Matrix4} m
     * @returns {Vector3}
     */
    transform(m) {
        const x = this.x, y = this.y, z = this.z;
        this.x = x*m.m00 + y*m.m10 + z*m.m20 + m.m30;
        this.y = x*m.m01 + y*m.m11 + z*m.m21 + m.m31;
        this.z = x*m.m02 + y*m.m12 + z*m.m22 + m.m32;
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    cross(other) {
        const x0 = this.x,  y0 = this.y,  z0 = this.z;
        const x1 = other.x, y1 = other.y, z1 = other.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {Vector3}
     */
    set(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    copy(other) {
        this.x = other.x;
        this.y = other.y;
        this.z = other.z;
        return this;
    }

    /**
     * @param {VectorObject3} other
     * @returns {Vector3}
     */
    neg_vector(other) {
        this.x = -other.x;
        this.y = -other.y;
        this.z = -other.z;
        return this;
    }

    /**
     * @param {VectorObject3} v0
     * @param {VectorObject3} v1
     * @returns {Vector3}
     */
    add_vectors(v0, v1) {
        this.x = v0.x+v1.x;
        this.y = v0.y+v1.y;
        this.z = v0.z+v1.z;
        return this;
    }

    /**
     * @param {VectorObject3} v0
     * @param {VectorObject3} v1
     * @returns {Vector3}
     */
    sub_vectors(v0, v1) {
        this.x = v0.x-v1.x;
        this.y = v0.y-v1.y;
        this.z = v0.z-v1.z;
        return this;
    }

    /**
     * @param {VectorObject3} v0
     * @param {VectorObject3} v1
     * @returns {Vector3}
     */
    mul_vectors(v0, v1) {
        this.x = v0.x*v1.x;
        this.y = v0.y*v1.y;
        this.z = v0.z*v1.z;
        return this;
    }

    /**
     * @param {VectorObject3} v0
     * @param {VectorObject3} v1
     * @returns {Vector3}
     */
    div_vectors(v0, v1) {
        this.x = v0.x/v1.x;
        this.y = v0.y/v1.y;
        this.z = v0.z/v1.z;
        return this;
    }

    /**
     * @param {VectorObject3} v0
     * @param {VectorObject3} v1
     * @param {number} t
     * @returns {Vector3}
     */
    lerp_vectors(v0, v1, t = 0) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;

        this.x = (x1 - x0) * t + x0;
        this.y = (y1 - y0) * t + y0;
        this.z = (z1 - z0) * t + z0;
        return this;
    }

    /**
     * @param {MatrixObject4} m
     * @param {VectorObject3} v
     * @returns {Vector3}
     */
    transform_vector(m, v) {
        const x = v.x, y = v.y, z = v.z;
        this.x = x*m.m00 + y*m.m10 + z*m.m20 + m.m30;
        this.y = x*m.m01 + y*m.m11 + z*m.m21 + m.m31;
        this.z = x*m.m02 + y*m.m12 + z*m.m22 + m.m32;
        return this;
    }

    /**
     * @param {VectorObject3} v0
     * @param {VectorObject3} v1
     * @returns {Vector3}
     */
    cross_vectors(v0, v1) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }
}

export class Vector4 {
    /**
     * @param {VectorObject4} v
     * @returns {Vector4}
     */
    static from(v)      { return new Vector4(v.x, v.y, v.z, v.w); }

    /**
     * @param {Number4} v
     * @returns {Vector4}
     */
    static fromArray(v) { return new Vector4(v[0], v[1], v[2], v[3]); }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @param {number} w
     * @constructor
     */
    constructor(x = 0, y = 0, z = 0, w = 1) {
        /**
         * @type {number}
         * @property x
         */
        this.x = x;

        /**
         * @type {number}
         * @property y
         */
        this.y = y;

        /**
         * @type {number}
         * @property z
         */
        this.z = z;

        /**
         * @type {number}
         * @property w
         */
        this.w = w;
    }

    *[Symbol.iterator]() {
        yield this.x;
        yield this.y;
        yield this.z;
        yield this.w;
    }

    /** @returns {Vector4} */
    clone() { return new Vector4(this.x, this.y, this.z, this.w); }

    /** @returns {VectorObject4} */
    toObject() { return { x: this.x, y: this.y, z: this.z, w: this.w }; }

    /** @returns {Number4} */
    toArray() { return [ this.x, this.y, this.z, this.w ]; }

    /**
     * @param {VectorObject4} other
     * @returns {number}
     */
    distance(other) {
        const x = other.x - this.x;
        const y = other.y - this.y;
        const z = other.z - this.z;
        return Math.sqrt(x*x + y*y + z*z);
    }

    /**
     * @returns {number}
     */
    magnitudeSquared() { return this.x*this.x + this.y*this.y + this.z*this.z; }

    /**
     * @returns {number}
     */
    magnitude() { return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z); }

    /**
     * @param {VectorObject4} other
     * @returns {number}
     */
    dot(other) { return this.x*other.x + this.y*other.y + this.z*other.z + this.w*other.w; }

    /**
     * @param {VectorObject4} other
     * @param {number} t
     * @returns {Vector4}
     */
    lerp(other, t = 0) {
        const x0 = this.x,  y0 = this.y,  z0 = this.z,  w0 = this.w;
        const x1 = other.x, y1 = other.y, z1 = other.z, w1 = other.w;

        return new Vector4((x1 - x0) * t + x0,
                           (y1 - y0) * t + y0,
                           (z1 - z0) * t + z0,
                           (w1 - w0) * t + w0);
    }

    /**
     * @returns {Vector4}
     */
    abs() {
        this.x = Math.abs(this.x);
        this.y = Math.abs(this.y);
        this.z = Math.abs(this.z);
        this.w = Math.abs(this.w);
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    min(other) {
        this.x = Math.min(this.x, other.x);
        this.y = Math.min(this.y, other.y);
        this.z = Math.min(this.z, other.z);
        this.w = Math.min(this.w, other.w);
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    max(other) {
        this.x = Math.max(this.x, other.x);
        this.y = Math.max(this.y, other.y);
        this.z = Math.max(this.z, other.z);
        this.w = Math.max(this.w, other.w);
        return this;
    }

    /**
     * @param {VectorObject4} v0
     * @param {VectorObject4} v1
     * @returns {Vector4}
     */
    clanp(v0, v1) {
        this.x = clamp(this.x, v0.x, v1.x);
        this.y = clamp(this.y, v0.y, v1.y);
        this.z = clamp(this.z, v0.z, v1.z);
        this.w = clamp(this.w, v0.w, v1.w);
        return this;
    }

    /**
     * @returns {Vector4}
     */
    neg() {
        this.x = -this.x;
        this.y = -this.y;
        this.z = -this.z;
        this.w = -this.w;
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    add(other) {
        this.x += other.x;
        this.y += other.y;
        this.z += other.z;
        this.w += other.w;
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    sub(other) {
        this.x -= other.x;
        this.y -= other.y;
        this.z -= other.z;
        this.w -= other.w;
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    mul(other) {
        this.x *= other.x;
        this.y *= other.y;
        this.z *= other.z;
        this.w *= other.w;
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    div(other) {
        this.x /= other.x;
        this.y /= other.y;
        this.z /= other.z;
        this.w /= other.w;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector4}
     */
    addScalar(scalar) {
        this.x += scalar;
        this.y += scalar;
        this.z += scalar;
        this.w += scalar;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector4}
     */
    subScalar(scalar) {
        this.x -= scalar;
        this.y -= scalar;
        this.z -= scalar;
        this.w -= scalar;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector4}
     */
    mulScalar(scalar) {
        this.x *= scalar;
        this.y *= scalar;
        this.z *= scalar;
        this.w *= scalar;
        return this;
    }

    /**
     * @param {number} scalar
     * @returns {Vector4}
     */
    divScalar(scalar) {
        this.x /= scalar;
        this.y /= scalar;
        this.z /= scalar;
        this.w /= scalar;
        return this;
    }

    /**
     * @param {number} length
     * @returns {Vector4}
     */
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

    /**
     * @returns {Vector4}
     */
    normalize() {
        const x0=this.x, y0=this.y, z0=this.z;
        const length_squared = x0*x0 + y0*y0 + z0*z0;

        if (length_squared == 0) {
            this.x = this.y = this.z = 0;
            return this;
        }

        const factor = 1 / Math.sqrt(length_squared);
        this.x = x0 * factor;
        this.y = y0 * factor;
        this.z = z0 * factor;
        return this;
    }

    /**
     * @param {MatrixObject4} matrix
     * @returns {Vector4}
     */
    transform(matrix) {
        const x = this.x, y = this.y, z = this.z, w = this.w;
        this.x = x*matrix.m00 + y*matrix.m10 + z*matrix.m20 + w*matrix.m30;
        this.y = x*matrix.m01 + y*matrix.m11 + z*matrix.m21 + w*matrix.m31;
        this.z = x*matrix.m02 + y*matrix.m12 + z*matrix.m22 + w*matrix.m32;
        this.w = x*matrix.m03 + y*matrix.m13 + z*matrix.m23 + w*matrix.m33;
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    cross(other) {
        const x0 = this.x,  y0 = this.y,  z0 = this.z;
        const x1 = other.x, y1 = other.y, z1 = other.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @param {number} w
     * @returns {Vector4}
     */
    set(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    copy(other) {
        this.x = other.x;
        this.y = other.y;
        this.z = other.z;
        this.w = other.w;
        return this;
    }

    /**
     * @param {VectorObject4} other
     * @returns {Vector4}
     */
    neg_vector(other) {
        this.x = -other.x;
        this.y = -other.y;
        this.z = -other.z;
        this.w = -other.w;
        return this;
    }

    /**
     * @param {VectorObject4} v0
     * @param {VectorObject4} v1
     * @returns {Vector4}
     */
    add_vectors(v0, v1) {
        this.x = v0.x+v1.x;
        this.y = v0.y+v1.y;
        this.z = v0.z+v1.z;
        this.w = v0.w+v1.w;
        return this;
    }

    /**
     * @param {VectorObject4} v0
     * @param {VectorObject4} v1
     * @returns {Vector4}
     */
    sub_vectors(v0, v1) {
        this.x = v0.x-v1.x;
        this.y = v0.y-v1.y;
        this.z = v0.z-v1.z;
        this.w = v0.w-v1.w;
        return this;
    }

    /**
     * @param {VectorObject4} v0
     * @param {VectorObject4} v1
     * @returns {Vector4}
     */
    mul_vectors(v0, v1) {
        this.x = v0.x*v1.x;
        this.y = v0.y*v1.y;
        this.z = v0.z*v1.z;
        this.w = v0.w*v1.w;
        return this;
    }

    /**
     * @param {VectorObject4} v0
     * @param {VectorObject4} v1
     * @returns {Vector4}
     */
    div_vectors(v0, v1) {
        this.x = v0.x/v1.x;
        this.y = v0.y/v1.y;
        this.z = v0.z/v1.z;
        this.w = v0.w/v1.w;
        return this;
    }

    /**
     * @param {VectorObject4} v0
     * @param {VectorObject4} v1
     * @param {number} t
     * @returns {Vector4}
     */
    lerp_vectors(v0, v1, t = 0) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z, w0 = v0.w;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z, w1 = v1.w;

        this.x = (x1 - x0) * t + x0;
        this.y = (y1 - y0) * t + y0;
        this.z = (z1 - z0) * t + z0;
        this.w = (w1 - w0) * t + w0;
        return this;
    }

    /**
     * @param {MatrixObject4} matrix
     * @param {VectorObject4} vector
     * @param {number} t
     * @returns {Vector4}
     */
    transform_vector(matrix, vector) {
        const x = vector.x, y = vector.y, z = vector.z, w = vector.w;
        this.x = x*matrix.m00 + y*matrix.m10 + z*matrix.m20 + w*matrix.m30;
        this.y = x*matrix.m01 + y*matrix.m11 + z*matrix.m21 + w*matrix.m31;
        this.z = x*matrix.m02 + y*matrix.m12 + z*matrix.m22 + w*matrix.m32;
        this.w = x*matrix.m03 + y*matrix.m13 + z*matrix.m23 + w*matrix.m33;
        return this;
    }

    /**
     * @param {VectorObject4} v0
     * @param {VectorObject4} v1
     * @param {number} t
     * @returns {Vector4}
     */
    cross_vectors(v0, v1) {
        const x0 = v0.x, y0 = v0.y, z0 = v0.z;
        const x1 = v1.x, y1 = v1.y, z1 = v1.z;
        this.x = y0 * z1 - y1 * z0;
        this.y = z0 * x1 - z1 * x0;
        this.z = x0 * y1 - x1 * y0;
        return this;
    }
}

export class Plane {
    /**
     * @param {Vector3} normal
     * @param {Vector3|Vector4} point
     * @returns {Plane}
     */
    static from_normal_and_point(normal, point) {
        const distance = normal.x*point.x + normal.y*point.y + normal.z*point.z;
        return new Plane(normal.x, normal.y, normal.z, distance);
    }

    /**
     * @param {Vector3|Vector4} point0
     * @param {Vector3|Vector4} point1
     * @param {Vector3|Vector4} point2
     * @returns {Plane}
     */
    static from_points(point0, point1, point2) {
        const x0      = point1.x - point0.x;
        const y0      = point1.y - point0.y;
        const z0      = point1.z - point0.z;
        const x1      = point2.x - point0.x;
        const y1      = point2.y - point0.y;
        const z1      = point2.z - point0.z;
        const cross_x = y0 * z1 - y1 * z0;
        const cross_y = z0 * x1 - z1 * x0;
        const cross_z = x0 * y1 - x1 * y0;
        const length  = Math.sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);

        if (length === 0) {
            return new Plane(0, 0, 0, 0);
        }

        const factor =  1 / length;
        const normal_x = cross_x * factor;
        const normal_y = cross_y * factor;
        const normal_z = cross_z * factor;
        const distance = normal_x*point0.x + normal_y*point0.y + normal_z*point0.z;

        return new Plane(normal_x, normal_y, normal_z, distance);
    }

    /**
     * @param {normal} normal_x
     * @param {normal} normal_y
     * @param {normal} normal_z
     * @param {normal} distance
     * @constructor
     */
    constructor(normal_x, normal_y, normal_z, distance) {
        this.normal_x = normal_x;
        this.normal_y = normal_y;
        this.normal_z = normal_z;
        this.distance = distance;
    }

    /**
     * @returns {Plane}
     */
    clone() {
        return new Plane(this.normal_x, this.normal_y, this.normal_z, this.distance);
    }

    /**
     * @param {Plane} other
     * @returns {Plane}
     */
    copy(other) {
        this.normal_x = other.normal_x;
        this.normal_y = other.normal_y;
        this.normal_z = other.normal_z;
        this.distance = other.distance;
        return this;
    }

    /**
     * @param {Plane} other
     * @returns {boolean}
     */
    equals(other) {
        const x = other.normal_x - this.normal_x;
        const y = other.normal_y - this.normal_y;
        const z = other.normal_z - this.normal_z;
        return Math.abs(this.distance - other.distance) <= EPSILON
        &&     Math.sqrt(x*x + y*y + z*z) <= EPSILON;
    }

    /**
     * @returns {Plane}
     */
    flip() {
        this.normal_x = -this.normal_x;
        this.normal_y = -this.normal_y;
        this.normal_z = -this.normal_z;
        this.distance = -this.distance;
        return this;
    }

    /**
     * @param {Vector3} point
     * @returns {number}
     */
    distance_to_point(point) {
        const dot_product = this.normal_x*point.x + this.normal_y*point.y + this.normal_z*point.z;
        return dot_product - this.distance;
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        return x * this.normal_x + y * this.normal_y + z * this.normal_z + this.distance;
    }

    /**
     * @param {Vector3} point
     * @returns {Vector3}
     */
    projection_of_point(point) {
        const a = this.distance_to_point(point);
        const x = point.x - a * this.normal_x;
        const y = point.y - a * this.normal_y;
        const z = point.z - a * this.normal_z;
        return new Vector3(x, y, z);
    }

    /**
     * @param {Matrix4} matrix
     * @param {Matrix4} inverse_matrix
     * @returns {Plane}
     */
    transform(matrix, inverse_matrix) {
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
}

export class Matrix4 {
    /**
     * @param {Matrix4} out
     * @returns {Matrix4}
     */
    static identity(out = new Matrix4()) {
        out.m00 = 1; out.m01 = 0; out.m02 = 0; out.m03 = 0;
        out.m10 = 0; out.m11 = 1; out.m12 = 0; out.m13 = 0;
        out.m20 = 0; out.m21 = 0; out.m22 = 1; out.m23 = 0;
        out.m30 = 0; out.m31 = 0; out.m32 = 0; out.m33 = 1;
        return out;
    }

    /**
     * @param {Vector3} v
     * @returns {Matrix4}
     */
    static scale(v) {
        const x = v.x, y = v.y, z = v.z;
        return new Matrix4(x, 0, 0, 0,
                           0, y, 0, 0,
                           0, 0, z, 0,
                           0, 0, 0, 1);
    }

    /**
     * @param {Vector3} v
     * @returns {Matrix4}
     */
    static translate(v) {
        const x = v.x, y = v.y, z = v.z;
        return new Matrix4(0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           x, y, z, 1);
    }

    /**
     * @param {number} angle
     * @returns {Matrix4}
     */
    static rotate_x(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(1, 0, 0, 0,
                           0, c,-s, 0,
                           0, s, c, 0,
                           0, 0, 0, 1);
    }

    /**
     * @param {number} angle
     * @returns {Matrix4}
     */
    static rotate_y(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(c, 0, s, 0,
                           0, 1, 0, 0,
                          -s, 0, c, 0,
                           0, 0, 0, 1);
    }

    /**
     * @param {number} angle
     * @returns {Matrix4}
     */
    static rotate_z(angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        return new Matrix4(c,-s, 0, 0,
                           s, c, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 1);
    }

    /**
     * @param {Vector3|Number3} axis
     * @param {number} angle
     * @returns {Matrix4}
     */
    static rotate(axis, angle) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        let x = axis.x, y = axis.y, z = axis.z;
        const l = Math.sqrt(x * x + y * y + z * z);
        const f = 1 / l;

        x *= f;
        y *= f;
        z *= f;

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

    /**
     * @param {number} m00
     * @param {number} m10
     * @param {number} m20
     * @param {number} m30
     * @param {number} m01
     * @param {number} m11
     * @param {number} m21
     * @param {number} m31
     * @param {number} m02
     * @param {number} m12
     * @param {number} m22
     * @param {number} m32
     * @param {number} m03
     * @param {number} m13
     * @param {number} m23
     * @param {number} m33
     * @constructor
     */
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

    /**
     * @returns {Matrix4}
     */
    clone() {
        return new Matrix4(this.m00, this.m10, this.m20, this.m30,
                           this.m01, this.m11, this.m21, this.m31,
                           this.m02, this.m12, this.m22, this.m32,
                           this.m03, this.m13, this.m23, this.m33);
    }

    /**
     * @param {MatrixObject4} other
     * @returns {Matrix4}
     */
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

    /**
     * @param {MatrixObject4} other
     * @returns {Matrix4}
     */
    transpose(other) {
        this.m00 = other.m00;
        this.m01 = other.m10;
        this.m02 = other.m20;
        this.m03 = other.m30;
        this.m10 = other.m01;
        this.m11 = other.m11;
        this.m12 = other.m21;
        this.m13 = other.m31;
        this.m20 = other.m02;
        this.m21 = other.m12;
        this.m22 = other.m22;
        this.m23 = other.m32;
        this.m30 = other.m03;
        this.m31 = other.m13;
        this.m32 = other.m23;
        this.m33 = other.m33;
        return this;
    }

    /**
     * @param {MatrixObject4} other
     * @returns {Matrix4}
     */
    mul(other) {
        return this.mul_matrices(this, other);
    }

    /**
     * @param {MatrixObject4} m0
     * @param {MatrixObject4} m1
     * @returns {Matrix4}
     */
    mul_matrices(m0, m1) {
        this.m00 = m0.m00 * m1.m00 + m0.m10 * m1.m01 + m0.m20 * m1.m02 + m0.m30 * m1.m03;
        this.m01 = m0.m00 * m1.m10 + m0.m10 * m1.m11 + m0.m20 * m1.m12 + m0.m30 * m1.m13;
        this.m02 = m0.m00 * m1.m20 + m0.m10 * m1.m21 + m0.m20 * m1.m22 + m0.m30 * m1.m23;
        this.m03 = m0.m00 * m1.m30 + m0.m10 * m1.m31 + m0.m20 * m1.m32 + m0.m30 * m1.m33;
        this.m10 = m0.m01 * m1.m00 + m0.m11 * m1.m01 + m0.m21 * m1.m02 + m0.m31 * m1.m03;
        this.m11 = m0.m01 * m1.m10 + m0.m11 * m1.m11 + m0.m21 * m1.m12 + m0.m31 * m1.m13;
        this.m12 = m0.m01 * m1.m20 + m0.m11 * m1.m21 + m0.m21 * m1.m22 + m0.m31 * m1.m23;
        this.m13 = m0.m01 * m1.m30 + m0.m11 * m1.m31 + m0.m21 * m1.m32 + m0.m31 * m1.m33;
        this.m20 = m0.m02 * m1.m00 + m0.m12 * m1.m01 + m0.m22 * m1.m02 + m0.m32 * m1.m03;
        this.m21 = m0.m02 * m1.m10 + m0.m12 * m1.m11 + m0.m22 * m1.m12 + m0.m32 * m1.m13;
        this.m22 = m0.m02 * m1.m20 + m0.m12 * m1.m21 + m0.m22 * m1.m22 + m0.m32 * m1.m23;
        this.m23 = m0.m02 * m1.m30 + m0.m12 * m1.m31 + m0.m22 * m1.m32 + m0.m32 * m1.m33;
        this.m30 = m0.m03 * m1.m00 + m0.m13 * m1.m01 + m0.m23 * m1.m02 + m0.m33 * m1.m03;
        this.m31 = m0.m03 * m1.m10 + m0.m13 * m1.m11 + m0.m23 * m1.m12 + m0.m33 * m1.m13;
        this.m32 = m0.m03 * m1.m20 + m0.m13 * m1.m21 + m0.m23 * m1.m22 + m0.m33 * m1.m23;
        this.m33 = m0.m03 * m1.m30 + m0.m13 * m1.m31 + m0.m23 * m1.m32 + m0.m33 * m1.m33;
        return this;
    }

    /**
     * @returns {number}
     */
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

    /**
     * @param {Matrix4} m
     * @returns {number}
     */
    transpose_inverse(m = this) {
        return this.inverse(m).transpose();
    }

    /**
     * @param {Matrix4} m
     * @returns {Matrix4}
     */
    inverse(m = this) {
        const subfactor0  = m.m22 * m.m33 - m.m23 * m.m32;
        const subfactor1  = m.m21 * m.m33 - m.m23 * m.m31;
        const subfactor2  = m.m21 * m.m32 - m.m22 * m.m31;
        const subfactor3  = m.m20 * m.m33 - m.m23 * m.m30;
        const subfactor4  = m.m20 * m.m32 - m.m22 * m.m30;
        const subfactor5  = m.m20 * m.m31 - m.m21 * m.m30;

        const subfactor6  = m.m02 * m.m13 - m.m03 * m.m12;
        const subfactor7  = m.m01 * m.m13 - m.m03 * m.m11;
        const subfactor8  = m.m01 * m.m12 - m.m02 * m.m11;
        const subfactor9  = m.m00 * m.m13 - m.m03 * m.m10;
        const subfactor10 = m.m00 * m.m12 - m.m02 * m.m10;
        const subfactor11 = m.m00 * m.m11 - m.m01 * m.m10;

        const cofactor00 = m.m11 * subfactor0  - m.m12 * subfactor1  + m.m13 * subfactor2;
        const cofactor01 = m.m12 * subfactor3  - m.m13 * subfactor4  - m.m10 * subfactor0;
        const cofactor02 = m.m10 * subfactor1  - m.m11 * subfactor3  + m.m13 * subfactor5;
        const cofactor03 = m.m11 * subfactor4  - m.m12 * subfactor5  - m.m10 * subfactor2;

        const cofactor10 = m.m02 * subfactor1  - m.m03 * subfactor2  - m.m01 * subfactor0;
        const cofactor11 = m.m00 * subfactor0  - m.m02 * subfactor3  + m.m03 * subfactor4;
        const cofactor12 = m.m01 * subfactor3  - m.m03 * subfactor5  - m.m00 * subfactor1;
        const cofactor13 = m.m00 * subfactor2  - m.m01 * subfactor4  + m.m02 * subfactor5;

        const cofactor20 = m.m31 * subfactor6  - m.m32 * subfactor7  + m.m33 * subfactor8;
        const cofactor21 = m.m32 * subfactor9  - m.m33 * subfactor10 - m.m30 * subfactor6;
        const cofactor22 = m.m30 * subfactor7  - m.m31 * subfactor9  + m.m33 * subfactor11;
        const cofactor23 = m.m31 * subfactor10 - m.m32 * subfactor11 - m.m30 * subfactor8;

        const cofactor30 = m.m22 * subfactor7  - m.m23 * subfactor8  - m.m21 * subfactor6;
        const cofactor31 = m.m20 * subfactor6  - m.m22 * subfactor9  + m.m23 * subfactor10;
        const cofactor32 = m.m21 * subfactor9  - m.m23 * subfactor11 - m.m20 * subfactor7;
        const cofactor33 = m.m20 * subfactor8  - m.m21 * subfactor10 + m.m22 * subfactor11;

        const determinant = m.m00 * cofactor00 + m.m01 * cofactor01 + m.m02 * cofactor02 + m.m03 * cofactor03;

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

    /**
     * @returns {boolean}
     */
    is_identity() {
        return this.m00 === 1 && this.m01 === 0 && this.m02 === 0 && this.m03 === 0
        &&     this.m10 === 0 && this.m11 === 1 && this.m12 === 0 && this.m13 === 0
        &&     this.m20 === 0 && this.m21 === 0 && this.m22 === 1 && this.m23 === 0
        &&     this.m30 === 0 && this.m31 === 0 && this.m32 === 0 && this.m33 === 1;
    }

    /**
     * @param {MatrixObject4} other
     * @returns {boolean}
     */
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

    /**
     * @param {VectorObject3} v
     * @param {Matrix4} m
     * @returns {Matrix4}
     */
    scale(v, m=this) {
        this.m00 = m.m00*v.x;
        this.m01 = m.m01*v.y;
        this.m02 = m.m02*v.z;
        this.m10 = m.m10*v.x;
        this.m11 = m.m11*v.y;
        this.m12 = m.m12*v.z;
        this.m20 = m.m20*v.x;
        this.m21 = m.m21*v.y;
        this.m22 = m.m22*v.z;
        this.m30 = m.m30*v.x;
        this.m31 = m.m31*v.y;
        this.m32 = m.m32*v.z;
        return this;
    }

    /**
     * @param {VectorObject3} v
     * @param {Matrix4} m
     * @returns {Matrix4}
     */
    translate(v, m=this) {
        this.m03 = m.m00*v.x + m.m01*v.y + m.m02*v.z + m.m03;
        this.m13 = m.m10*v.x + m.m11*v.y + m.m12*v.z + m.m13;
        this.m23 = m.m20*v.x + m.m21*v.y + m.m22*v.z + m.m23;
        this.m33 = m.m30*v.x + m.m31*v.y + m.m32*v.z + m.m33;
        return this;
    }

    /**
     * @param {number} angle
     * @param {Matrix4} m
     * @returns {Matrix4}
     */
    rotate_x(angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m01 = m.m01*c - m.m02*s;
        this.m02 = m.m02*c + m.m01*s;
        this.m11 = m.m11*c - m.m12*s;
        this.m12 = m.m12*c + m.m11*s;
        this.m21 = m.m21*c - m.m22*s;
        this.m22 = m.m22*c + m.m21*s;
        this.m31 = m.m31*c - m.m32*s;
        this.m32 = m.m32*c + m.m31*s;
        return this;
    }

    /**
     * @param {number} angle
     * @param {Matrix4} m
     * @returns {Matrix4}
     */
    rotate_y(angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = m.m00*c + m.m02*s;
        this.m02 = m.m02*c - m.m00*s;
        this.m10 = m.m10*c + m.m12*s;
        this.m12 = m.m12*c - m.m10*s;
        this.m20 = m.m20*c + m.m22*s;
        this.m22 = m.m22*c - m.m20*s;
        this.m30 = m.m30*c + m.m32*s;
        this.m32 = m.m32*c - m.m30*s;
        return this;
    }

    /**
     * @param {number} angle
     * @param {Matrix4} m
     * @returns {Matrix4}
     */
    rotate_z(angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);

        this.m00 = m.m00*c - m.m01*s;
        this.m01 = m.m01*c + m.m00*s;
        this.m10 = m.m10*c - m.m11*s;
        this.m11 = m.m11*c + m.m10*s;
        this.m20 = m.m20*c - m.m21*s;
        this.m21 = m.m21*c + m.m20*s;
        this.m30 = m.m30*c - m.m31*s;
        this.m31 = m.m31*c + m.m30*s;
        return this;
    }

    /**
     * @param {VectorObject3} axis
     * @param {number} angle
     * @param {Matrix4} m
     * @returns {Matrix4}
     */
    rotate(axis, angle, m=this) {
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        const one_minus_c = 1 - c;

        let x = axis.x, y = axis.y, z = axis.z;
        const l = Math.sqrt(x * x + y * y + z * z);
        const f = 1 / l;

        x *= f;
        y *= f;
        z *= f;

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

        this.m00 = m.m10 * b10 + m.m20 * b20 + m.m00 * b00; // m00
        this.m01 = m.m00 * b01 + m.m20 * b21 + m.m10 * b11; // m01
        this.m02 = m.m00 * b02 + m.m10 * b12 + m.m20 * b22; // m02
        this.m10 = m.m11 * b10 + m.m21 * b20 + m.m01 * b00; // m10
        this.m11 = m.m01 * b01 + m.m21 * b21 + m.m11 * b11; // m11
        this.m12 = m.m01 * b02 + m.m11 * b12 + m.m21 * b22; // m12
        this.m20 = m.m12 * b10 + m.m22 * b20 + m.m02 * b00; // m20
        this.m21 = m.m02 * b01 + m.m22 * b21 + m.m12 * b11; // m21
        this.m22 = m.m02 * b02 + m.m12 * b12 + m.m22 * b22; // m22
        this.m30 = m.m13 * b10 + m.m23 * b20 + m.m03 * b00; // m30
        this.m31 = m.m03 * b01 + m.m23 * b21 + m.m13 * b11; // m31
        this.m32 = m.m03 * b02 + m.m13 * b12 + m.m23 * b22; // m32
        return this;
    }
}

export class BoundingBox {
    /**
     * @param {SignedDistanceField[]} fields
     * @returns {BoundingBox}
     */
    static enclose_distance_fields(fields) {
        const bounding_box = fields[0].bounding_box.clone();

        for (let idx = 1; idx < fields.length; ++idx) {
            bounding_box.enclose_box(fields[idx].bounding_box);
        }

        return bounding_box;
    }

    /**
     * @param {SignedDistanceField[]} fields
     * @returns {BoundingBox}
     */
    static intersecting_distance_fields(fields) {
        const bounding_box = fields[0].bounding_box.clone();

        for (let idx = 1; idx < fields.length; ++idx) {
            bounding_box.intersect_box(fields[idx].bounding_box);
        }

        return bounding_box;
    }

    /**
     * @param {BoundingBox} box0
     * @param {BoundingBox} box1
     * @returns {BoundingBox}
     */
    static create_enclosing_box(box0, box1) {
        return new BoundingBox(
            Math.min(box0.min_x, box1.min_x),
            Math.min(box0.min_y, box1.min_y),
            Math.min(box0.min_z, box1.min_z),
            Math.max(box0.max_x, box1.max_x),
            Math.max(box0.max_y, box1.max_y),
            Math.max(box0.max_z, box1.max_z),
        );
    }

    /**
     * @param {BoundingBox} box0
     * @param {BoundingBox} box1
     * @returns {BoundingBox}
     */
    static create_intersecting_box(box0, box1) {
        return new BoundingBox(
            Math.max(box0.min_x, box1.min_x),
            Math.max(box0.min_y, box1.min_y),
            Math.max(box0.min_z, box1.min_z),
            Math.min(box0.max_x, box1.max_x),
            Math.min(box0.max_y, box1.max_y),
            Math.min(box0.max_z, box1.max_z),
        );
    }

    /**
     * @param {number} min_x
     * @param {number} min_y
     * @param {number} min_z
     * @param {number} max_x
     * @param {number} max_y
     * @param {number} max_z
     * @constructor
     */
    constructor(min_x = 0, min_y = 0, min_z = 0, max_x = 0, max_y = 0, max_z = 0) {
        /** @type {number} */
        this.min_x = min_x;

        /** @type {number} */
        this.min_y = min_y;
        
        /** @type {number} */
        this.min_z = min_z;
        
        /** @type {number} */
        this.max_x = max_x;
        
        /** @type {number} */
        this.max_y = max_y;
        
        /** @type {number} */
        this.max_z = max_z;
    }

    /**
     * @returns {number}
     */
    size_x() { return this.max_x - this.min_x; }

    /**
     * @returns {number}
     */
    size_y() { return this.max_y - this.min_y; }

    /**
     * @returns {number}
     */
    size_z() { return this.max_z - this.min_z; }

    /**
     * @returns {number}
     */
    largest_side() { return Math.max(this.size_x(), this.size_y(), this.size_z()); }

    /**
     * @returns {number}
     */
    smallest_side() { return Math.min(this.size_x(), this.size_y(), this.size_z()); }

    /**
     * @returns {number}
     */
    smallest_grid_size() {
        const size_x = this.max_x - this.min_x;
        const size_y = this.max_y - this.min_y;
        const size_z = this.max_z - this.min_z;
        return greatest_common_factor(size_x,
               greatest_common_factor(size_y, size_z));
    }

    /**
     * @param {number} grid_size
     * @returns {number}
     */
    minimal_grid_size(grid_size) {
        const size_x = this.max_x - this.min_x;
        const size_y = this.max_y - this.min_y;
        const size_z = this.max_z - this.min_z;

        if (size_x == size_y && size_x == size_z) {
            const outside_volume = size_x / grid_size;
            const inside_volume = Math.floor(outside_volume);
            return (inside_volume / outside_volume) * grid_size;
        }

        const smallest_grid_size = greatest_common_factor(
            greatest_common_factor(size_x, size_y), size_z);

        if (grid_size < smallest_grid_size) {
            return Math.min(Math.floor(smallest_grid_size / grid_size), 1) * grid_size;
        }

        return Math.max(Math.floor(grid_size / smallest_grid_size), 1) * smallest_grid_size;
    }

    /**
     * @param {number} step_size
     * @returns {BoundingBox}
     */
    optimize_smallest_grid_size(step_size) {
        const center_x = (this.max_x + this.min_x) / 2;
        const center_y = (this.max_y + this.min_y) / 2;
        const center_z = (this.max_z + this.min_z) / 2;
        const size_x = step_size * Math.round((this.max_x - this.min_x) / step_size);
        const size_y = step_size * Math.round((this.max_y - this.min_y) / step_size);
        const size_z = step_size * Math.round((this.max_z - this.min_z) / step_size);

        this.min_x = center_x - size_x / 2;
        this.min_y = center_y - size_y / 2;
        this.min_z = center_z - size_z / 2;
        this.max_x = center_x + size_x / 2;
        this.max_y = center_y + size_y / 2;
        this.max_z = center_z + size_z / 2;
        return this;
    }

    /**
     * @param {Matrix4} matrix
     * @returns {BoundingBox}
     */
    transform(matrix) {
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

    /**
     * @returns {BoundingBox}
     */
    clone() {
        return new BoundingBox(
            this.min_x,
            this.min_y,
            this.min_z,
            this.max_x,
            this.max_y,
            this.max_z,
        );
    }

    /**
     * @param {number} size
     * @returns {BoundingBox}
     */
    grow(size) {
        this.min_x -= size;
        this.min_y -= size;
        this.min_z -= size;
        this.max_x += size;
        this.max_y += size;
        this.max_z += size;
        return this;
    }

    /**
     * @param {number} factor
     * @returns {BoundingBox}
     */
    scale(factor) {
        this.min_x *= factor;
        this.min_y *= factor;
        this.min_z *= factor;
        this.max_x *= factor;
        this.max_y *= factor;
        this.max_z *= factor;
        return this;
    }

    /**
     * @param {BoundingBox} other
     * @returns {BoundingBox}
     */
    enclose_box(other) {
        this.min_x = Math.min(this.min_x, other.min_x);
        this.min_y = Math.min(this.min_y, other.min_y);
        this.min_z = Math.min(this.min_z, other.min_z);
        this.max_x = Math.max(this.max_x, other.max_x);
        this.max_y = Math.max(this.max_y, other.max_y);
        this.max_z = Math.max(this.max_z, other.max_z);
        return this;
    }

    /**
     * @param {BoundingBox} other
     * @returns {BoundingBox}
     */
    intersect_box(other) {
        this.min_x = Math.max(this.min_x, other.min_x);
        this.min_y = Math.max(this.min_y, other.min_y);
        this.min_z = Math.max(this.min_z, other.min_z);
        this.max_x = Math.min(this.max_x, other.max_x);
        this.max_y = Math.min(this.max_y, other.max_y);
        this.max_z = Math.min(this.max_z, other.max_z);
        return this;
    }

    /**
     * @param {Polygon[]} polygons
     * @returns {Polygon[]}
     */
    to_polygons(polygons = []) {
        const position0 = new Vector3(this.min_x, this.min_y, this.min_z);
        const position1 = new Vector3(this.max_x, this.min_y, this.min_z);
        const position2 = new Vector3(this.max_x, this.max_y, this.min_z);
        const position3 = new Vector3(this.min_x, this.max_y, this.min_z);
        const position4 = new Vector3(this.min_x, this.min_y, this.max_z);
        const position5 = new Vector3(this.max_x, this.min_y, this.max_z);
        const position6 = new Vector3(this.max_x, this.max_y, this.max_z);
        const position7 = new Vector3(this.min_x, this.max_y, this.max_z);
        polygons.push(new Polygon([ position3, position2, position1, position0 ]));
        polygons.push(new Polygon([ position4, position5, position6, position7 ]));
        polygons.push(new Polygon([ position0, position1, position5, position4 ]));
        polygons.push(new Polygon([ position2, position3, position7, position6 ]));
        polygons.push(new Polygon([ position1, position2, position6, position5 ]));
        polygons.push(new Polygon([ position4, position7, position3, position0 ]));
        return polygons;
    }
}

export class Vertex {
    /**
     * @param {Vector3} position
     * @param {Vector3} normal
     * @param {Vector4} colour
     * @param {Vector3} texture_coords
     * @param {number} texture_id
     * @constructor
     */
    constructor(
        position = new Vector3(0, 0, 0),
        normal = new Vector3(0, 0, 0),
        colour = new Vector4(0, 0, 0, 0),
        texture_coords = new Vector3(0, 0, 0),
        texture_id = 0,
    ) {
        /** @type {Vector3} */
        this.position = position;
        /** @type {Vector3} */
        this.normal = normal;
        /** @type {Vector4} */
        this.colour = colour;
        /** @type {Vector3} */
        this.texture_coords = texture_coords;
        /** @type {number} */
        this.texture_id = texture_id;
    }
}

export class PolygonFace {
    /**
     * @param {number[]} indices
     * @param {Plane} plane
     * @constructor
     */
    constructor(indices, plane) {
        this.indices = indices;
        this.plane = plane;
    }
}

export class PolygonMesh {
    /**
     * @param {Vertex[]} vertices
     * @param {PolygonFace[]} faces
     * @constructor
     */
    constructor(vertices, faces) {
        this.vertices = vertices;
        this.faces = faces;
    }

    /**
     * @returns {string}
     */
    serialize_to_obj() {
        const string = "";

        for (let idx = 0; idx < this.vertices.length; ++idx) {
            const position = this.vertices[idx].position;
            string += `v\t${position.x}\t${position.y}\t${position.z}\n`;
        }

        for (let idx = 0; idx < this.faces.length; ++idx) {
            const indices = this.faces[idx].indices;
            string += `v\t${indices.join('\t')}\n`;
        }

        return string;
    }
    
    /**
     * @returns {Uint8Array}
     */
    serialize_to_stl() {
        const buffer = new ArrayBuffer(84 + 50 * this.faces.length);
        const view = new DataView(buffer);

        view.setUint32(80, this.faces.length, true);

        for (let offset = 84, idx = 0; idx < this.faces.length; ++idx) {
            const face = this.faces[idx];
            
            if (face.indices.length > 3) {
                return null;
            }

            const vertex0 = this.vertices[face.indices[0]];
            const vertex1 = this.vertices[face.indices[1]];
            const vertex2 = this.vertices[face.indices[2]];

            view.setFloat32(offset, face.plane.normal_x, true); offset += 4;
            view.setFloat32(offset, face.plane.normal_y, true); offset += 4;
            view.setFloat32(offset, face.plane.normal_z, true); offset += 4;

            view.setFloat32(offset, vertex0.position.x, true); offset += 4;
            view.setFloat32(offset, vertex0.position.y, true); offset += 4;
            view.setFloat32(offset, vertex0.position.z, true); offset += 4;

            view.setFloat32(offset, vertex1.position.x, true); offset += 4;
            view.setFloat32(offset, vertex1.position.y, true); offset += 4;
            view.setFloat32(offset, vertex1.position.z, true); offset += 4;

            view.setFloat32(offset, vertex2.position.x, true); offset += 4;
            view.setFloat32(offset, vertex2.position.y, true); offset += 4;
            view.setFloat32(offset, vertex2.position.z, true); offset += 4;

            view.setUint16(offset, 0, true); offset += 2;
        }

        return new Uint8Array(buffer);
    }
}

export class Polygon {
    constructor(positions) {
        this.positions = positions;
    }
}

export class SignedDistanceField {
    static temp_vector0 = new Vector3();

    /**
     * @param {BoundingBox} bounding_box
     * @constructor
     */
    constructor(bounding_box = new BoundingBox()) {
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_UNKNOWN;
        /** @type {BoundingBox} */
        this.bounding_box = bounding_box;
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        throw new Error("calculate_signed_distance is unimplemented!");
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @param {Vector3} out_vector
     * @returns {Vector3}
     */
    calculate_gradient(x, y, z, out_vector = new Vector3(0, 0, 0)) {
        const distance_max_x = this.calculate_signed_distance(x + EPSILON, y, z);
        const distance_min_x = this.calculate_signed_distance(x - EPSILON, y, z);
        const distance_max_y = this.calculate_signed_distance(x, y + EPSILON, z);
        const distance_min_y = this.calculate_signed_distance(x, y - EPSILON, z);
        const distance_max_z = this.calculate_signed_distance(x, y, z + EPSILON);
        const distance_min_z = this.calculate_signed_distance(x, y, z - EPSILON);
        out_vector.set((distance_max_x - distance_min_x) * INV_TWO_EPSILON,
                       (distance_max_y - distance_min_y) * INV_TWO_EPSILON,
                       (distance_max_z - distance_min_z) * INV_TWO_EPSILON);
        return out_vector;
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @param {Vector3} out_vector
     * @returns {Vector3}
     */
    calculate_closest_surface_vector(x, y, z, out_vector = new Vector3(0, 0, 0)) {
        const distance = this.calculate_signed_distance(x, y, z);
        const delta    = this.calculate_gradient(x, y, z, out_vector).setLength(distance);
        return out_vector.set(x - delta.x, y - delta.y, z - delta.z);
    }

    /**
     * @param {SignedDistanceField[]} others
     * @returns {SignedDistanceFieldUnion}
     */
    union(... others) {
        return new SignedDistanceFieldUnion(this, ... others);
    }

    /**
     * @param {SignedDistanceField[]} others
     * @returns {SignedDistanceFieldDifference}
     */
    difference(... others) {
        return new SignedDistanceFieldDifference(this, ... others);
    }

    /**
     * @param {SignedDistanceField[]} others
     * @returns {SignedDistanceFieldIntersection}
     */
    intersection(... others) {
        return new SignedDistanceFieldIntersection(this, ... others);
    }

    /**
     * @param {number} smoothness
     * @param {SignedDistanceField[]} others
     * @returns {SignedDistanceFieldUnionSmooth}
     */
    union_smooth(smoothness, ... others) {
        return new SignedDistanceFieldUnionSmooth(smoothness, this, ... others);
    }

    /**
     * @param {number} smoothness
     * @param {SignedDistanceField[]} others
     * @returns {SignedDistanceFieldDifferenceSmooth}
     */
    difference_smooth(smoothness, ... others) {
        return new SignedDistanceFieldDifferenceSmooth(smoothness, this, ... others);
    }

    /**
     * @param {number} smoothness
     * @param {SignedDistanceField[]} others
     * @returns {SignedDistanceFieldIntersectionSmooth}
     */
    intersection_smooth(smoothness, ... others) {
        return new SignedDistanceFieldIntersectionSmooth(smoothness, this, ... others);
    }

    /**
     * @param {Matrix4} matrix
     * @returns {SignedDistanceFieldTransform}
     */
    transform(matrix) {
        return new SignedDistanceFieldTransform(this, matrix);
    }

    /**
     * @param {number} angle
     * @returns {SignedDistanceFieldTransform}
     */
    rotate_x(angle) {
        return this.transform(Matrix4.rotate_x(angle));
    }

    /**
     * @param {number} angle
     * @returns {SignedDistanceFieldTransform}
     */
    rotate_y(angle) {
        return this.transform(Matrix4.rotate_y(angle));
    }

    /**
     * @param {number} angle
     * @returns {SignedDistanceFieldTransform}
     */
    rotate_z(angle) {
        return this.transform(Matrix4.rotate_z(angle));
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {SignedDistanceFieldTranslate}
     */
    translate(x, y, z) {
        return new SignedDistanceFieldTranslate(x, y, z, this);
    }

    /**
     * @param {number} radius
     * @returns {SignedDistanceFieldOffset}
     */
    offset(radius) {
        return new SignedDistanceFieldOffset(this, radius);
    }

    /**
     * @param {number} step_count
     * @returns {PolygonMesh}
     */
    surface_nets(step_count = 200) {
        const step_size             = this.bounding_box.largest_side() / step_count;
        const bounding_box          = this.bounding_box.clone().optimize_smallest_grid_size(step_size);
        const cube_size             = bounding_box.minimal_grid_size(step_size);
        const sample_count_x        = Math.round(bounding_box.size_x() / cube_size) + 2;
        const sample_count_y        = Math.round(bounding_box.size_y() / cube_size) + 2;
        const sample_count_z        = Math.round(bounding_box.size_z() / cube_size) + 2;
        const half_size             = cube_size / 2;

        console.log(cube_size);
        console.log(this.bounding_box.largest_side() / cube_size);

        let index = 0;
        const vertices = [], faces = [];
        for (let idx0 = 0, x = bounding_box.min_x - half_size; idx0 < sample_count_x; ++idx0, x += cube_size) {
            for (let idx1 = 0, y = bounding_box.min_y - half_size; idx1 < sample_count_y; ++idx1, y += cube_size) {
                for (let idx2 = 0, z = bounding_box.min_z - half_size; idx2 < sample_count_z; ++idx2, z += cube_size, ++index) {
                    const sample_origin = this.calculate_signed_distance(x, y, z);
                    const sample_x      = this.calculate_signed_distance(x + cube_size, y, z);
                    const sample_y      = this.calculate_signed_distance(x, y + cube_size, z);
                    const sample_z      = this.calculate_signed_distance(x, y, z + cube_size);

                    if ((sample_origin < 0) != (sample_x < 0)) {
                        const index = vertices.length;
                        vertices.push(new Vertex(new Vector3(x + half_size, y - half_size, z - half_size)));
                        vertices.push(new Vertex(new Vector3(x + half_size, y + half_size, z - half_size)));
                        vertices.push(new Vertex(new Vector3(x + half_size, y + half_size, z + half_size)));
                        vertices.push(new Vertex(new Vector3(x + half_size, y - half_size, z + half_size)));
                        if (sample_origin < 0) {
                            faces.push(new PolygonFace([ index + 0, index + 1, index + 2 ], null));
                            faces.push(new PolygonFace([ index + 2, index + 3, index + 0 ], null));
                        }
                        else {
                            faces.push(new PolygonFace([ index + 2, index + 1, index + 0 ], null));
                            faces.push(new PolygonFace([ index + 0, index + 3, index + 2 ], null));
                        }
                    }

                    if ((sample_origin < 0) != (sample_y < 0)) {
                        const index = vertices.length;
                        vertices.push(new Vertex(new Vector3(x - half_size, y + half_size, z + half_size)));
                        vertices.push(new Vertex(new Vector3(x + half_size, y + half_size, z + half_size)));
                        vertices.push(new Vertex(new Vector3(x + half_size, y + half_size, z - half_size)));
                        vertices.push(new Vertex(new Vector3(x - half_size, y + half_size, z - half_size)));
                        if (sample_origin < 0) {
                            faces.push(new PolygonFace([ index + 0, index + 1, index + 2 ], null));
                            faces.push(new PolygonFace([ index + 2, index + 3, index + 0 ], null));
                        }
                        else {
                            faces.push(new PolygonFace([ index + 2, index + 1, index + 0 ], null));
                            faces.push(new PolygonFace([ index + 0, index + 3, index + 2 ], null));
                        }
                    }

                    if ((sample_origin < 0) != (sample_z < 0)) {
                        const index = vertices.length;
                        vertices.push(new Vertex(new Vector3(x - half_size, y - half_size, z + half_size)));
                        vertices.push(new Vertex(new Vector3(x + half_size, y - half_size, z + half_size)));
                        vertices.push(new Vertex(new Vector3(x + half_size, y + half_size, z + half_size)));
                        vertices.push(new Vertex(new Vector3(x - half_size, y + half_size, z + half_size)));
                        if (sample_origin < 0) {
                            faces.push(new PolygonFace([ index + 0, index + 1, index + 2 ], null));
                            faces.push(new PolygonFace([ index + 2, index + 3, index + 0 ], null));
                        }
                        else {
                            faces.push(new PolygonFace([ index + 2, index + 1, index + 0 ], null));
                            faces.push(new PolygonFace([ index + 0, index + 3, index + 2 ], null));
                        }
                    }
                }
            }
        }

        for (let idx = 0; idx < vertices.length; ++idx) {
            const position = vertices[idx].position;
            this.calculate_closest_surface_vector(position.x, position.y, position.z, position);
            this.calculate_closest_surface_vector(position.x, position.y, position.z, position);
            this.calculate_closest_surface_vector(position.x, position.y, position.z, position);
            this.calculate_closest_surface_vector(position.x, position.y, position.z, position);
        }

        for (let idx = 0; idx < faces.length; ++idx) {
            const vertex0 = vertices[faces[idx].indices[0]];
            const vertex1 = vertices[faces[idx].indices[1]];
            const vertex2 = vertices[faces[idx].indices[2]];
            vertex0.normal = this.calculate_gradient(vertex0.position.x, vertex0.position.y, vertex0.position.z).normalize();
            vertex1.normal = this.calculate_gradient(vertex1.position.x, vertex1.position.y, vertex1.position.z).normalize();
            vertex2.normal = this.calculate_gradient(vertex2.position.x, vertex2.position.y, vertex2.position.z).normalize();
            faces[idx].plane = Plane.from_points(vertex0.position, vertex1.position, vertex2.position);
        }

        return new PolygonMesh(vertices, faces);
    }
}

export class SignedDistanceFieldSphere extends SignedDistanceField {
    constructor(center_x = 0, center_y = 0, center_z = 0, radius = 5) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_SPHERE;
        /** @type {number} */
        this.center_x = center_x;
        /** @type {number} */
        this.center_y = center_y;
        /** @type {number} */
        this.center_z = center_z;
        /** @type {number} */
        this.radius   = radius;
        /** @type {BoundingBox} */
        this.bounding_box = new BoundingBox(
            this.center_x - radius,
            this.center_y - radius,
            this.center_z - radius,
            this.center_x + radius,
            this.center_y + radius,
            this.center_z + radius,
        );
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        const dx = x - this.center_x;
        const dy = y - this.center_y;
        const dz = z - this.center_z;
        return Math.sqrt(dx*dx + dy*dy + dz*dz) - this.radius;
    }
}

export class SignedDistanceFieldCylinder extends SignedDistanceField {
    constructor(center_x = 0, center_y = 0, center_z = 0, height = 20, radius = 10) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_CYLINDER;
        /** @type {number} */
        this.center_x = center_x;
        /** @type {number} */
        this.center_y = center_y;
        /** @type {number} */
        this.center_z = center_z;
        /** @type {number} */
        this.radius      = radius;
        /** @type {number} */
        this.half_height = height / 2;
        /** @type {BoundingBox} */
        this.bounding_box = new BoundingBox(
            this.center_x - this.radius,
            this.center_y - this.radius,
            this.center_z - this.half_height,
            this.center_x + this.radius,
            this.center_y + this.radius,
            this.center_z + this.half_height,
        );
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        const dx = x - this.center_x;
        const dy = y - this.center_y;
        const dz = z - this.center_z;
        const d0 = Math.abs(Math.sqrt(dx*dx + dy*dy)) - this.radius;
        const d1 = Math.abs(dz) - this.half_height;
        const a0 = Math.max(d0,0.0);
        const a1 = Math.max(d1,0.0);
        return Math.min(Math.max(d0,d1),0.0) + Math.sqrt(a0*a0 + a1*a1);
    }
}

export class SignedDistanceFieldPlane extends SignedDistanceField {
    constructor(plane) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_PLANE;
        /** @type {Plane} */
        this.plane = plane;
        /** @type {BoundingBox} */
        this.bounding_box = new BoundingBox(0, 0, 0, 0, 0, 0);
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        return this.plane.calculate_signed_distance(x, y, z);
    }
}

export class SignedDistanceFieldCube extends SignedDistanceField {
    constructor(center_x = 0, center_y = 0, center_z = 0, size = 10, corner_radius = 0) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_CUBE;
        /** @type {number} */
        this.center_x = center_x;
        /** @type {number} */
        this.center_y = center_y;
        /** @type {number} */
        this.center_z = center_z;
        /** @type {number} */
        this.half_size = size / 2;
        /** @type {number} */
        this.corner_radius = corner_radius;
        /** @type {BoundingBox} */
        this.bounding_box = new BoundingBox(
            this.center_x - this.half_size,
            this.center_y - this.half_size,
            this.center_z - this.half_size,
            this.center_x + this.half_size,
            this.center_y + this.half_size,
            this.center_z + this.half_size
        );
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        const pos_x = Math.abs(x - this.center_x) - this.half_size + this.corner_radius;
        const pos_y = Math.abs(y - this.center_y) - this.half_size + this.corner_radius;
        const pos_z = Math.abs(z - this.center_z) - this.half_size + this.corner_radius;
        const side = Math.min(Math.max(pos_x, pos_y, pos_z), 0.0);
        const dx = Math.max(pos_x, 0);
        const dy = Math.max(pos_y, 0);
        const dz = Math.max(pos_z, 0);
        return Math.sqrt(dx*dx + dy*dy + dz*dz) + side - this.corner_radius;
    }
}

export class SignedDistanceFieldBox extends SignedDistanceField {
    constructor(center_x = 0, center_y = 0, center_z = 0, size_x = 20, size_y = 20, size_z = 20, corner_radius = 0) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_BOX;
        /** @type {number} */
        this.center_x = center_x;
        /** @type {number} */
        this.center_y = center_y;
        /** @type {number} */
        this.center_z = center_z;
        /** @type {number} */
        this.half_size_x = 0.5 * size_x;
        /** @type {number} */
        this.half_size_y = 0.5 * size_y;
        /** @type {number} */
        this.half_size_z = 0.5 * size_z;
        /** @type {number} */
        this.corner_radius = corner_radius;
        /** @type {BoundingBox} */
        this.bounding_box  = new BoundingBox(
            this.center_x - this.half_size_x,
            this.center_y - this.half_size_y,
            this.center_z - this.half_size_z,
            this.center_x + this.half_size_x,
            this.center_y + this.half_size_y,
            this.center_z + this.half_size_z,
        );
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        const pos_x = Math.abs(x - this.center_x) - this.half_size_x + this.corner_radius;
        const pos_y = Math.abs(y - this.center_y) - this.half_size_y + this.corner_radius;
        const pos_z = Math.abs(z - this.center_z) - this.half_size_z + this.corner_radius;
        const side = Math.min(Math.max(pos_x, pos_y, pos_z), 0.0);
        const dx = Math.max(pos_x, 0);
        const dy = Math.max(pos_y, 0);
        const dz = Math.max(pos_z, 0);
        return Math.sqrt(dx*dx + dy*dy + dz*dz) + side - this.corner_radius;
    }
}

export class SignedDistanceFieldUnion extends SignedDistanceField {
    /**
     * @param  {SignedDistanceField[]} objects 
     * @constructor
     */
    constructor(... objects) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_UNION;
        /** @type  {SignedDistanceField[]} */
        this.objects = objects;
        /** @type {BoundingBox} */
        this.bounding_box = BoundingBox.enclose_distance_fields(this.objects);
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        let distance = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            distance = Math.min(distance, this.objects[idx].calculate_signed_distance(x, y, z));
        }

        return distance;
    }
}

export class SignedDistanceFieldDifference extends SignedDistanceField {
    /**
     * @param  {SignedDistanceField[]} objects 
     * @constructor
     */
    constructor(... objects) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_DIFFERENCE;
        /** @type  {SignedDistanceField[]} */
        this.objects = objects;
        /** @type {BoundingBox} */
        this.bounding_box = this.objects[0].bounding_box.clone();
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        let distance = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            distance = Math.max(distance, -this.objects[idx].calculate_signed_distance(x, y, z));
        }

        return distance;
    }
}

export class SignedDistanceFieldIntersection extends SignedDistanceField {
    /**
     * @param  {SignedDistanceField[]} objects 
     * @constructor
     */
    constructor(... objects) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_INTERSECTION;
        /** @type  {SignedDistanceField[]} */
        this.objects = objects;
        /** @type {BoundingBox} */
        this.bounding_box = BoundingBox.intersecting_distance_fields(this.objects);
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        let distance = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            distance = Math.max(distance, this.objects[idx].calculate_signed_distance(x, y, z));
        }

        return distance;
    }
}

export class SignedDistanceFieldUnionSmooth extends SignedDistanceField {
    constructor(smoothness, ... objects) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_UNION_SMOOTH;
        this.smoothness = smoothness;
        this.objects = objects;
        /** @type {BoundingBox} */
        this.bounding_box = BoundingBox.enclose_distance_fields(this.objects);
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        let distance0 = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            let distance1 = this.objects[idx].calculate_signed_distance(x, y, z);
            let h = clamp(0.5 + 0.5 * (distance1 - distance0) / this.smoothness, 0, 1);
            let m = distance1 + (distance0 - distance1) * h;
            distance0 = m - this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldDifferenceSmooth extends SignedDistanceField {
    constructor(smoothness, ... objects) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_DIFFERENCE_SMOOTH;
        this.smoothness = smoothness;
        this.objects = objects;
        /** @type {BoundingBox} */
        this.bounding_box = this.objects[0].bounding_box.clone();
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        let distance0 = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            let distance1 = this.objects[idx].calculate_signed_distance(x, y, z);
            let h = clamp(0.5 - 0.5 * (distance1 + distance0) / this.smoothness, 0, 1);
            let m = distance0 + (-distance1 - distance0) * h;
            distance0 = m + this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldIntersectionSmooth extends SignedDistanceField {
    constructor(smoothness, ... objects) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_INTERSECTION_SMOOTH;
        this.smoothness = smoothness;
        this.objects = objects;
        /** @type {BoundingBox} */
        this.bounding_box = BoundingBox.intersecting_distance_fields(this.objects);
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        let distance0 = this.objects[0].calculate_signed_distance(x, y, z);

        for (let idx = 1; idx < this.objects.length; ++idx) {
            let distance1 = this.objects[idx].calculate_signed_distance(x, y, z);
            let h = clamp(0.5 - 0.5 * (distance1 - distance0) / this.smoothness, 0, 1);
            let m = distance1 + (distance0 - distance1) * h;
            distance0 = m + this.smoothness * h * (1 - h);
        }

        return distance0;
    }
}

export class SignedDistanceFieldTranslate extends SignedDistanceField {
    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @param {SignedDistanceField} field
     * @constructor
     */
    constructor(x, y, z, field) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_TRANSLATE;
        /** @type {number} */
        this.translate_x = x;
        /** @type {number} */
        this.translate_y = y;
        /** @type {number} */
        this.translate_z = z;
        /** @type {SignedDistanceField} */
        this.field       = field;
        /** @type {BoundingBox} */
        this.bounding_box = new BoundingBox(
            field.bounding_box.min_x + x,
            field.bounding_box.min_y + y,
            field.bounding_box.min_z + z,
            field.bounding_box.max_x + x,
            field.bounding_box.max_y + y,
            field.bounding_box.max_z + z,
        );
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        return this.field.calculate_signed_distance(x - this.translate_x,
                                                    y - this.translate_y,
                                                    z - this.translate_z);
    }
}

export class SignedDistanceFieldTransform extends SignedDistanceField {
    /** @type {Vector3} */
    static temp_vector = new Vector3(0, 0, 0);

    /**
     * @param {SignedDistanceField} field
     * @param {Matrix4} matrix 
     * @constructor
     */
    constructor(field, matrix) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_TRANSFORM;
        /** @type {SignedDistanceField} */
        this.field = field;
        /** @type {Matrix4} */
        this.matrix = matrix.clone();
        /** @type {Matrix4} */
        this.inverse_matrix = matrix.inverse();
        /** @type {BoundingBox} */
        this.bounding_box = this.field.bounding_box.clone().transform(this.matrix);
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        SignedDistanceFieldTransform.temp_vector.set(x, y, z);
        SignedDistanceFieldTransform.temp_vector.transform(this.inverse_matrix);
        return this.field.calculate_signed_distance(
            SignedDistanceFieldTransform.temp_vector.x,
            SignedDistanceFieldTransform.temp_vector.y,
            SignedDistanceFieldTransform.temp_vector.z);
    }
}

export class SignedDistanceFieldOffset extends SignedDistanceField {
    constructor(object, offset) {
        super();
        /** @type {number} */
        this.kind = SIGNED_DISTANCE_FIELD_KIND_OFFSET;
        /** @type {number} */
        this.radius = offset;
        /** @type {SignedDistanceField} */
        this.field = object;
        /** @type {BoundingBox} */
        this.bounding_box = new BoundingBox(
            this.field.bounding_box.min_x - this.radius,
            this.field.bounding_box.min_y - this.radius,
            this.field.bounding_box.min_z - this.radius,
            this.field.bounding_box.max_x + this.radius,
            this.field.bounding_box.max_y + this.radius,
            this.field.bounding_box.max_z + this.radius,
        );
    }

    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @returns {number}
     */
    calculate_signed_distance(x, y, z) {
        return this.field.calculate_signed_distance(x, y, z) - this.radius;
    }
}

export function cube({ center = [0, 0, 0], size = 10, corner_radius = 0 }) {
    return new SignedDistanceFieldCube(center[0], center[1], center[2], size, corner_radius);
}

export function box({ center = [0, 0, 0], size = [ 20, 20, 20 ], corner_radius = 0 }) {
    return new SignedDistanceFieldBox(center[0], center[1], center[2], size[0], size[1], size[2], corner_radius);
}

export function sphere({ center = [0, 0, 0], radius = 5 }) {
    return new SignedDistanceFieldSphere(center[0], center[1], center[2], radius);
}

export function cylinder({ center = [ 0, 0, 0 ], height = 20, radius = 10 }) {
    return new SignedDistanceFieldCylinder(center[0], center[1], center[2], height, radius);
}

export function plane({ center = [0, 0, 0] }) {
    return new SignedDistanceFieldPlane(center[0], center[1], center[2], height, radius);
}

export function torus({ center = [0, 0, 0] }) {
    return new SignedDistanceFieldPlane(center[0], center[1], center[2], height, radius);
}
