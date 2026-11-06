# Ask-a-question assistant for the report site (planning)

Started 2026-11-04, after a senior team member asked for an LLM on the site: a user asks a question, the site
answers with supporting plots. Status (2026-11-06): **the AI-service team answered** (summary below); the design is
updated to match. Nothing built yet; next steps need decisions (see "Decisions needed").

## Design (agreed in discussion)

- The model **routes and summarizes; it does not analyze.** It turns a question into settings for the site's
  existing, checked tools (cost map explorer, peak tool, region, platforms, data viewer), runs them, and writes a
  short answer from the numbers they return, with the actual views embedded and a link that reproduces each one.
- The browser keeps doing the data work; only the question and a few summary numbers go to the model.
- A small server-side piece is needed to hold the API credentials (GitHub Pages can't).
- Audience: **the team** first; possibly collaborating oceanographers at other institutions later (same expertise).
  Not public-facing.
- Likely questions: the ten core analyses' themes (the boss's original interests), e.g. global mean T/S against
  depth, where and at what depth the cost is high, split vs unsplit, cost against time, per-bin maps, zonal means,
  outliers, and attribution of peaks.

## Built regardless (useful with or without an LLM; also the assistant's "actions") — done 2026-11-04

1. `questions.html` ("What do you want to know?"): 21 questions from the ten analyses' themes, each a preset view
   (`diagnostics/profile_analyses/questions.py`; the list is the natural starting set of assistant actions).
2. Every view as a link: the report (`?src=&ct=&an=&band=&part=`), the cost map explorer (`?src=&ct=&V=&part=
   &depth=&years=&range=&same=&region=&log=&max=`, plus `peak=` to run the peak tool) and the data viewer
   (`?src=&ref=&V=&k=|dmode=&g=&y0=&y1=&platform=&view=`). An assistant can answer with these links and embed them.
3. Simple / expert mode in the explorer.

## Message for the AI-service contact (fill in the brackets)

**Subject:** Using [the internal AI service] from a small team web tool: what's possible?

Hi [name],

Our team has built a web tool for exploring ocean profile data (temperature and salinity) and how well it fits
our ocean model. The data are open, and it's currently a static site. People find the many options hard to
navigate, so we'd like to let a user type a question, e.g. "Where does the peak in Argo cost near 2400 m come
from?", and have the page answer with the relevant plots.

The model's job would be narrow: turn the question into settings for our existing analysis tools, then
summarize the numbers those tools return. The analysis itself stays in our code; the model wouldn't see raw data,
only the question and a few summary numbers. Users would be our team, possibly a few collaborating oceanographers
later; it's not public-facing.

Could you help with these questions?
1. **Access:** can a web page or small server we run call [the service] through an API? If so, how do we get
   credentials, and are there usage limits or costs to the team?
2. **Function / tool calling:** do the available models support calling functions we define (structured outputs
   or tool use), or only plain chat?
3. **Where it can run:** where are we allowed to host the small server piece that holds the credentials?
   Something like an internal VM, a cloud function in our team's AWS account, or a service you provide? Our site
   is on GitHub Pages now; can a page there call it, or does the whole tool need to be hosted internally?
4. **Users outside [the institution]:** if we later share the tool with collaborators at other institutions,
   can they use it, or is it internal-only?
5. **Data policy:** the questions and numbers are about open, published ocean data. Is there anything we should
   know about what may or may not be sent, logged, or retained?
6. **Models and docs:** which models are available, and is there documentation or an example of a team
   integrating it into their own tool?

Happy to show a quick demo of what we have. Thanks!

### Why these questions
- **Tool calling** decides whether the model can drive the tools reliably.
- **Where it can run** likely decides the architecture; "internal only" may mean moving the whole site behind the
  institution's login (which would also remove the secure-network blocking seen with R2).
- Saying the data are open and only summaries are sent eases the policy conversation.

## Answers from the AI-service team (2026-11-06, summarized)

| Question | Answer | What it means for us |
|---|---|---|
| 1. Access, cost | Two services: the **GenAI Enterprise Managed API** (Enterprise Tech; models through Amazon Bedrock; open to anyone; subscribe at goto.jpl.nasa.gov/GenAI-API-Service) and **Onyx** (CAE; closed to new users; new instance from 2026-11-16 for existing users and two priority projects). Managed API subscriptions: **Individual** (your JPL login; for prototyping) or **Application** (a JPL application account; for a site others use). $11/month plus tokens, on a project charge number (may need the CAM's approval) or, Individual only, the timecard. A monthly cap ($100 by default) stops calls when reached. Help: #peer-genai-managed-api, genai-api-admin@jpl.nasa.gov. | Use the **Managed API**: Individual to prototype, Application for the team. Cost is small at our scale. Need a charge number. |
| 2. Tool calling | Both support it. The Managed API speaks the OpenAI chat format and handles tool calling routinely. | Our design (the model picks among our actions) is a normal pattern there. |
| 3. Where it runs | The credential must live on a small server we run, never in the page (policy). Page → our proxy (holds the credential, limits requests per user) → model. Both services are on the JPL network, so the proxy should be too (e.g. a JPL VM). Application tokens last 12 hours, no refresh: the proxy logs in again. The Managed API examples include a **Streamlit** app (the whole tool in one internal Python app). **A public github.com Pages site can't reach an internal proxy; github.jpl.nasa.gov/pages sites are internal and fine.** Don't open a proxy to the internet. A production app may need approval (the service runs under JPL IT Plan 95; apps under another security plan need that plan owner's and the CISO office's sign-off; our information system owner knows which plan applies). | **Our site is on public github.com**, so as it stands the page can't call an internal proxy. Either move (a copy of) the site to github.jpl.nasa.gov Pages, or run the whole assistant as an internal app. |
| 4. Outside users | Neither service lets people outside JPL sign in. Foreign nationals at JPL can subscribe by email request. Whether a JPL tool may serve model answers to outside collaborators: ask genai-api-admin first. | **Team-only (JPL) for now.** External collaborators: a separate policy question, later. |
| 5. Data | Managed API: GovCloud (default) cleared up to export-controlled; Commercial cleared up to CUI, with more models and slightly lower prices. Neither allows personal information. For public data, Commercial is reasonable. Retention of prompts and answers isn't documented: ask. The page **must** show "This content was generated by AI and may contain mistakes. Validate for accuracy." and, if the model answers without our data or sources, say so. Throttle on our side; under load expect 429 / 503 errors. Send the model the chart's current state and the question, not the dataset; keep visitors' names out. | Matches our design: the model gets the question and the tools' summary numbers, never data. Add the required notice; add throttling and retry on 429/503. Ask about retention. |
| 6. Models, docs | Managed API (GovCloud, October): Claude Opus 5.5, Sonnet 5.5, Haiku 5.5, plus Amazon models (Commercial has more). Onyx: Claude Sonnet 4.5, 5, 5.5, Haiku 5.5. **Start with Haiku 5.5**; move to Sonnet 5.5 if it struggles. Docs: genai.jpl.nasa.gov/gen-ai-tools/managed-service (user guide, acceptable use policy); examples github.jpl.nasa.gov/FutureTechGenAI/genai_api; JPL GenAI catalog fiddle.jpl.nasa.gov:8081/catalog. | Haiku 5.5 for routing questions to actions; Sonnet 5.5 if summaries need more care. |

**Also a policy rule:** when AI enhances an existing feature, the original way must stay. The assistant sits **next
to** the buttons, not instead of them, which is already our plan.

## Design, updated

```
team member's browser (JPL network)
  page on github.jpl.nasa.gov Pages  ──question + current view──▶  proxy on a JPL VM  ──▶  Managed API (Haiku 5.5)
        ▲    runs the chosen actions in the browser                  holds the credential,          (tool calling)
        │    (the existing checked tools)                            per-user limits, logs in
        └──── answer text + links/views ◀── summary numbers ◀────── every 12 h, retries 429/503
```
1. The user asks a question in a box next to the existing controls.
2. The page sends the proxy the question and the current view (class, cost type, variable, depth, years), never data.
3. The model chooses actions from a fixed list (built from questions.html and the view links): e.g. open_view,
   run_peak(range), set_region(name or box), list_platforms, open_profile(platform).
4. The page runs them with the existing tools, and sends the model back only their summary numbers (as the peak
   panel prints them).
5. The model writes a short answer from those numbers; the page shows it with the views and links, and the
   required AI notice. Every number in the answer is one the tools produced, and the links reproduce them.

Alternative (simpler to host, more to rebuild): the whole site as an internal Streamlit app.

## Decisions needed
1. **Hosting.** Mirror the site to **github.jpl.nasa.gov Pages** (keeps the public github.com site as it is) and add
   a **JPL VM** for the proxy, or build an internal Streamlit app. Recommended: the mirror + VM (reuses everything).
2. **Charge number** for a Managed API subscription (and the CAM's approval if needed); **Individual** to prototype,
   **Application** once the team uses it.
3. **GovCloud or Commercial.** The data are public, so Commercial is allowed and slightly cheaper; GovCloud is the
   default. Either works for Haiku / Sonnet 5.5.
4. **Who to ask about production approval** (the project's information system owner), before team-wide use.

## Questions to send back (short)
- How long are prompts and answers retained, and where (for both GovCloud and Commercial)?
- May a tool on github.jpl.nasa.gov Pages call a proxy on a JPL VM? (the answer suggests yes)
- For a team-internal tool: is the Individual subscription acceptable for a small prototype others try?

## Next
Once the decisions above are made: a prototype with an Individual subscription, a local proxy, and the actions
for the ten analyses' questions; then the internal mirror and the VM.
