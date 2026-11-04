# Ask-a-question assistant for the report site (planning)

Started 2026-11-04, after a senior team member asked for an LLM on the site: a user asks a question, the site
answers with supporting plots. Status: **waiting on the institution's AI-service team** (message below).

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

## Next
When the answers come back: a concrete design (where the server runs, which model, the list of actions with their
inputs and outputs), then a prototype that answers the ten analyses' questions through the actions.
