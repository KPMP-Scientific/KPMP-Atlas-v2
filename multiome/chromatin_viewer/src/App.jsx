/* eslint-disable */
import * as React from 'react';
import {
  ThemeProvider,
  StylesProvider,
  createGenerateClassName,
} from '@material-ui/core';
import Viewer from './Viewer';

const generateClassName = createGenerateClassName({
	disableGlobal: false, // Class names need to be deterministic,
});

function App() {
	return (
		<StylesProvider generateClassName={generateClassName}>
			<ThemeProvider>
				<Viewer />
			</ThemeProvider>
		</StylesProvider>
	);
}

export default App;