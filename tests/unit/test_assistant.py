"""The assistant's question loop with a stand-in model: tool calls are run, their links collected, and the answer
carries the required AI notice. No network."""
from ecco_insitu.diagnostics.assistant import server, tools


class FakeModel:
    model = 'fake'

    def __init__(self):
        self.calls = []

    def chat(self, messages, schemas):
        self.calls.append(list(messages))       # a copy: the loop keeps adding to the list
        if len(self.calls) == 1:     # first: ask for one tool
            return {'choices': [{'message': {'role': 'assistant', 'content': None, 'tool_calls': [
                {'id': 't1', 'type': 'function', 'function': {'name': 'list_contents', 'arguments': '{}'}}]}}]}
        return {'choices': [{'message': {'role': 'assistant', 'content': 'There are 2 classes.'}}], 'usage': {}}


class FakeAnalysis:
    sources, cost_types = ['PFL', 'XBT_WOD'], ['vs_clim']


def test_answer_runs_tools_and_adds_notice(monkeypatch):
    seen = []
    monkeypatch.setattr(tools, 'call', lambda A, name, args: (seen.append(name), ({'obs_classes': A.sources},
                        [{'page': 'report', 'params': {'an': '02'}, 'caption': 'x'}]))[1])
    m = FakeModel()
    out = server.answer(FakeAnalysis(), m, 'what classes are there?', {'page': 'report', 'src': 'PFL'})
    assert seen == ['list_contents']
    assert out['answer'] == 'There are 2 classes.'
    assert out['notice'] == server.AI_NOTICE
    assert out['views'] == [{'page': 'report', 'params': {'an': '02'}, 'caption': 'x'}]
    assert m.calls[1][-1]['role'] == 'tool'                       # the tool's result went back to the model
    assert '"src": "PFL"' in m.calls[0][0]['content']            # the page's view is in the instructions


def test_unknown_tool_is_reported_not_raised():
    res, views = tools.call(None, 'no_such_tool', {})
    assert 'error' in res and views == []
