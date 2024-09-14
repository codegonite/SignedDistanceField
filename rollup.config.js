import terser from '@rollup/plugin-terser';

export default {
	input: 'src/index.js',
	output: [
        { file: 'build/SignedDistanceField.cjs',     format: 'cjs' },
        { file: 'build/SignedDistanceField.js',      format: 'esm' },
        { file: 'build/SignedDistanceField.min.cjs', format: 'cjs', plugins: [terser()] },
        { file: 'build/SignedDistanceField.min.js',  format: 'esm', plugins: [terser()] },
    ]
};