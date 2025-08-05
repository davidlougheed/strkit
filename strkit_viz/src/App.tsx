import { Link, Outlet } from "react-router-dom";

const App = () => {
  return (
    <div id="strkit">
      <header>
        <h1>STRkit Browser</h1>
        <nav>
          <ul>
            <li><Link to="/">Overview</Link></li>
            <li><Link to="/loci">Loci</Link></li>
          </ul>
        </nav>
      </header>
      <Outlet />
    </div>
  );
};

export default App;
