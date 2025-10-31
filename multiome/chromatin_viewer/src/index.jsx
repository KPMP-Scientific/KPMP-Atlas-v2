import React from 'react';
import ReactDOM from 'react-dom';
import { NuqsAdapter } from 'nuqs/adapters/react';
import App from './App';
import './index.css';

ReactDOM.render(
	<NuqsAdapter>
		<App />
	</NuqsAdapter>,
	document.getElementById('root')
);
