class Interface:

    def __init__(self, conn):
        self.canvas = conn.ui.stock_canvas
        self.panel = self.canvas.add_panel()
        self.panel.rect_transform.size = (380, 100)
        self.panel.rect_transform.position = (-560, -475)
        self.phase_current = self.panel.add_text('Phase Current: ')
        self.phase_next = self.panel.add_text('Phase Next: ')
        self.phase_current_time = self.panel.add_text('(T-{:0>5d})'.format(0))
        self.next_action = self.panel.add_text('Next Action: N/A')
        self.message = self.panel.add_text('Message: ')
        self.button = self.panel.add_button("Execute Action")

        self.phase_current.rect_transform.position = (-80, 30)
        self.phase_next.rect_transform.position = (-80, 10)
        self.phase_current_time.rect_transform.position = (65, 30)
        self.next_action.rect_transform.position = (-50, -10)
        self.message.rect_transform.position = (-50, -30)

        self.phase_current.rect_transform.size = (200, 20)
        self.phase_next.rect_transform.size = (200, 20)
        self.phase_current_time.rect_transform.size = (80, 20)
        self.next_action.rect_transform.size = (260, 20)
        self.message.rect_transform.size = (260, 20)
        self.button.rect_transform.size = (90, 90)
        self.button.rect_transform.position = (140, 0)

        self.button_clicked = conn.add_stream(getattr, self.button, 'clicked')

    def update_phase_current(self, s):
        self.phase_current.content = 'Current: ' + s

    def update_phase_next(self, s):
        self.phase_next.content = 'Next: ' + s

    def update_phase_current_time(self, s):
        self.phase_current_time.content = '(T-{:0>5d})'.format(s)

    def update_next_action(self, s):
        self.next_action.content = 'Action: ' + s

    def update_message(self, s):
        self.message.content = 'Message: ' + s
