import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  base: "/KPMP-Atlas-v2/chromatin-viewer/",
  build: {
    target: "esnext",
  },
  plugins: [
    react({
      jsxRuntime: "classic",
    }),
  ],
  // To enable .js files that contain JSX to be imported.
  // Reference: https://github.com/vitest-dev/vitest/issues/1564
  esbuild: {
    loader: "tsx",
    include: /src\/.*\.[tj]sx?$/,
    exclude: [],
  },
});
