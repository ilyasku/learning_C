project('First Try Harmonic Oscillator', 'c')
cc = meson.get_compiler('c')
mathdep = cc.find_library('m')
executable('harmonic', ['src/main.c', 'src/harmonic_oscillator.c'],
	  dependencies: mathdep)

checkdep = dependency('check')

test_harmonic_executable = executable('test_executable', ['tests/test_harmonic_oscillator.c',
							  'src/harmonic_oscillator.c'],
				     dependencies: [checkdep, mathdep])
test('test harmonic oscillator', test_harmonic_executable)
