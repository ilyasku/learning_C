project('First Try Harmonic Oscillator', 'c')
executable('harmonic', ['src/main.c', 'src/harmonic_oscillator.c'])

checkdep = dependency('check')

test_harmonic_executable = executable('test_executable', ['tests/test_harmonic_oscillator.c',
							  'src/harmonic_oscillator.c'],
				     dependencies: checkdep)
test('test harmonic oscillator', test_harmonic_executable)
