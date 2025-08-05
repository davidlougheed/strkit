import { createRoot } from "react-dom/client";
import { createBrowserRouter, Navigate, RouterProvider } from "react-router-dom";

import App from "./App";
import Overview from "./Overview";
import Locus from "./Locus";

const router = createBrowserRouter([
  {
    path: "/",
    Component: App,
    children: [
      {
        index: true,
        Component: Overview,
      },
      {
        path: "loci",
        Component: () => <Navigate to="/loci/0" />,
      },
      {
        path: "loci/:locusIdx",
        Component: Locus,
      },
    ],
  },
]);

document.addEventListener("DOMContentLoaded", () => {
  const root = createRoot(document.getElementById("root"));

  root.render(
      <RouterProvider router={router} />
  );
});
