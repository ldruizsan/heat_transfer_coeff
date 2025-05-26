# Engineering Heat Transfer Chatbot - Development Roadmap

This document outlines the planned development phases for the Engineering Heat Transfer Chatbot. The goal is to create a tool that assists undergraduate engineering students in identifying appropriate empirical correlations, understanding their applicability, and calculating heat transfer coefficients.

## Core Philosophy

*   **Educational Value:** Prioritize features that help students learn *why* a correlation is used, not just get a number.
*   **Modularity:** Maintain a clear separation between the calculation engine, NLU, dialogue management, and UI.
*   **Iterative Development:** Start with core functionality and incrementally add features and complexity.

## Phase 0: Foundation - Calculation Engine Enhancement (Largely Complete)

*   **[✓] Task:** Refine `app_with_fn_handler.py` (the calculation engine).
    *   Nusselt handler functions (`_nusselt_..._detailed`) return a rich dictionary including:
        *   `Nu_value`
        *   `correlation_name`
        *   `correlation_formula_str`
        *   `conditions_applied` (Re, Pr, flow regime, etc.)
        *   `explanation_notes` (justification, assumptions)
        *   `warnings`
    *   Main `calculate_heat_transfer_coefficient` function returns `h_value`, `Nu_details` dictionary, and `error` message.
*   **[ ] Task:** Develop comprehensive unit tests for `app_with_fn_handler.py`.
    *   Cover various geometries, thermal conditions, and flow regimes.
    *   Test edge cases and Pr/Re boundary conditions for correlations.
    *   Verify the accuracy of the returned explanatory details.

## Phase 1: Core Chatbot - CLI, Basic NLU, and Fluid Properties

*Goal: A command-line chatbot that can answer a fully specified query and provide detailed explanations.*

*   **[ ] Task: Command-Line Interface (CLI) - `chatbot_ui.py` / `main_chatbot.py`**
    *   Develop the main application loop for user interaction.
    *   Handle basic input and output.
*   **[ ] Task: Simple Natural Language Understanding (NLU) - `nlu_parser.py`**
    *   Implement keyword spotting and regular expressions to extract:
        *   Geometry type (e.g., "flat plate", "cylinder")
        *   Thermal condition (e.g., "isothermal", "constant heat flux")
        *   Key parameters (e.g., velocity, characteristic length)
        *   Fluid type and temperature (e.g., "air at 300K")
    *   Convert extracted entities into the structured format required by the calculation engine.
*   **[ ] Task: Fluid Property Database - `fluid_database.py`**
    *   Create an initial database (e.g., Python dictionary, JSON) for common fluids (Air, Water) at a few standard temperatures.
    *   Function to look up `rho`, `mu`, `cp`, `k` based on fluid name and approximate temperature.
*   **[ ] Task: Response Generation - `response_generator.py`**
    *   Develop functions to take the detailed dictionary from `calculate_heat_transfer_coefficient` and format it into a clear, explanatory textual response.
    *   Include all relevant details: chosen correlation, formula, justification, calculated values, and warnings.
*   **[ ] Task: Integration**
    *   Connect all Phase 1 components: UI -> NLU -> Fluid DB -> Calculation Engine -> Response Generator -> UI.
    *   Handle cases where the user provides all necessary information in a single query.

## Phase 2: Enhanced NLU & Basic Dialogue Management

*Goal: A more robust chatbot that can ask clarifying questions if the user's query is incomplete.*

*   **[ ] Task: Improve NLU Parser (`nlu_parser.py`)**
    *   Handle more variations in user phrasing.
    *   Consider integrating a lightweight NLP library (e.g., `spaCy` for Named Entity Recognition) if regex becomes too complex.
*   **[ ] Task: Basic Dialogue Management (in `chatbot_ui.py`)**
    *   Implement "slot filling": If essential parameters are missing from the initial query (e.g., fluid velocity, characteristic length), the chatbot should prompt the user for them.
    *   Maintain a simple state to track required vs. provided information for the current query.
*   **[ ] Task: Expand Fluid Database (`fluid_database.py`)**
    *   Add more fluids and a wider range of temperatures.
    *   Consider simple interpolation for temperatures not directly listed.
*   **[ ] Task: Refine Response Generation (`response_generator.py`)**
    *   Handle conversational prompts for missing information.
    *   Improve error messages for NLU failures or missing fluid properties.

## Phase 3: Advanced Features & Usability

*Goal: A more comprehensive and user-friendly engineering assistant.*

*   **[ ] Task: Unit Handling & Conversion**
    *   Allow users to input values with common engineering units (e.g., length in cm/mm/in, temperature in °C/°F).
    *   Implement internal conversion to SI units for calculations.
*   **[ ] Task: Contextual Help & Examples**
    *   Implement commands like "help", "supported geometries", "example query for flat plate".
*   **[ ] Task: Expand Geometries & Correlations (in `app_with_fn_handler.py` and NLU)**
    *   Incrementally add more scenarios:
        *   Sphere in crossflow
        *   Internal pipe flow (laminar, turbulent)
        *   Natural convection (vertical plate, horizontal cylinder)
    *   Update NLU to recognize new geometries and related parameters.
    *   Update response generator for new scenarios.
*   **[ ] Task: Enhanced Error Handling & Feedback**
    *   Provide more specific feedback from the calculation engine (e.g., "Reynolds number [value] is too low for the selected turbulent correlation, which typically applies for Re > [limit]").
*   **[ ] Task: Input Summary & Confirmation**
    *   Before performing a calculation, have the chatbot summarize the interpreted parameters and ask for user confirmation, especially for complex queries.
*   **[ ] Task: Logging**
    *   Implement logging for user queries (anonymized if necessary) and chatbot responses for debugging and identifying common issues or desired features.

## Phase 4: Potential Future Enhancements (Beyond Core Scope)

*These are ideas for longer-term development if the core chatbot is successful and widely used.*

*   **[ ] Task: Web Interface**
    *   Develop a simple web UI (e.g., using Flask, Django, or Streamlit) for broader accessibility.
*   **[ ] Task: Graphical Output**
    *   Consider generating simple diagrams of the geometry being analyzed.
    *   Plot relevant dimensionless numbers or show property variations if temperature-dependent properties are implemented.
*   **[ ] Task: Advanced NLU/Dialogue**
    *   Explore more sophisticated NLU/dialogue management frameworks like Rasa or cloud-based AI services for more natural and complex conversations.
*   **[ ] Task: User Accounts & History (for web version)**
    *   Allow users to save common fluid properties, previous queries, or preferences.
*   **[ ] Task: Integration with External Fluid Property Databases**
    *   Connect to more comprehensive databases (e.g., CoolProp) for a wider range of fluids and more accurate property calculations.
*   **[ ] Task: Support for Mixed Boundary Layer Conditions**
    *   Implement correlations for flat plates where flow starts laminar and transitions to turbulent.

## Versioning

*   **v0.1:** Completion of Phase 0 (Solid Calculation Engine).
*   **v0.5:** Completion of Phase 1 (Basic functional CLI chatbot).
*   **v1.0:** Completion of Phase 2 (CLI chatbot with dialogue for parameter elicitation).
*   **v1.x - v2.0:** Incremental releases incorporating features from Phase 3.

---
This roadmap is a living document and may be updated as the project progresses and new requirements or ideas emerge.
```

This roadmap provides a structured plan. Each phase builds upon the previous one, allowing for iterative development and testing. Remember to adapt it as you go based on your progress and any new insights you gain!